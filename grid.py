#!/usr/bin/env python

import tkinter as tk

import cv2  # resize: si integral to match rorig
import matplotlib.pyplot as plt  # for color scale
import nibabel as nib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # numpy plots
from matplotlib.figure import Figure
from PIL import Image, ImageTk

from siarray import Scout, SIArray
import lcmodel


# ## create new circle+box methods for tk.Canvas
def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x - r, y - r, x + r, y + r, **kwargs)


def _create_box(self, x, y, r, **kwargs):
    return self.create_rectange(x - r, y - r, x + r, y + r, **kwargs)


tk.Canvas.create_circle = _create_circle
tk.Canvas.create_box = _create_box


def npimg(a):
    return ImageTk.PhotoImage(image=Image.fromarray(a))


class ROI:
    """store roi"""

    def __init__(self, label: str, xy):
        self.roi = label
        self.xy = xy[0:2]

    def update(self, event):
        """update on coords on mouse click
        N.B. works for axial. need to reorient for cor or sag"""
        self.xy = [event.x, event.y]

    def move(self, x=0, y=0, dim=None):
        "adjust x and/or y position (arrow press). wrap around dimensions"
        if dim is None:
            dim = [100000, 100000]

        self.xy = [(self.xy[0] + x) % dim[0], (self.xy[1] + y) % dim[1]]

    def sid3(self, res_edge):
        "return (y, x) reverse direction"
        # [int(res_edge - p) for p in reverse(self.xy)]
        return (int(self.xy[1]), int(res_edge - self.xy[0]))

    def label(self, res_edge):
        pos = self.sid3(res_edge)
        return f"{self.roi} {pos[0]} {pos[1]}"


class App(tk.Frame):
    def __init__(self, master=None, roixy_list=None, ref=None, si=None, sires=24):
        super().__init__(master)
        self.master = master
        self.master.title("MRSI Coord Placer")

        # INIT
        self.fnames = {"t1": ref, "si": si, "gm": None}
        self.imgs = {"ax-": None, "ax0": None, "ax+": None, "si": None}
        self.i_curroi = tk.IntVar(self)
        self.coords = None  # [ROI('roi1',[x, y])...]
        # update later with integral of siarray
        self.siIntg = None
        self.scout = None  # holds scout resolution

        # DIMS
        # TODO: these are for 7T MRSI. could be setting somewhere
        #       instead of hardcoded
        self.pixdim = (216, 216, 99)
        self.voxdim = (9, 9, 10)
        self.sires = (sires, sires)

        # VIS
        self.pack()
        self.create_widgets()

        # set by set_coords during load
        # needed for label update and to match box outline to listbox bg
        self.hexcolors = []
        self.set_coords(roixy_list)

        # autoload when start if t1 and si are defined
        self.load()

    def inc_roi_selected(self, step=1):
        "set roi number: global var and pos in listbox"
        n = self.roiselect.size()
        next_roi = (self.i_curroi.get() + step) % n
        self.i_curroi.set(next_roi)
        self.roiselect.selection_clear(0, n)
        self.roiselect.selection_set(next_roi)
        self.roiselect.see(next_roi)
        # need to update colors. redraw images and boxes
        self.draw_images()
        self.add_coords()

    def move_roi(self, x=0, y=0, step=1):
        "move the selected roi by one grid element"
        roi = self.i_curroi.get()
        self.coords[roi].move(x * step, y * step)
        # print(f"moving {roi} {x},{y}: {self.coords[roi].xy}")
        self.update()

    def create_widgets(self):

        # TODO: these could be top bar menu item?
        self.btnfrm = tk.Frame(highlightbackground="blue", highlightthickness=2)
        self.loadbtn = tk.Button(self.btnfrm, text="load", command=self.load)
        self.loadbtn.pack(side="left")
        self.savebtn = tk.Button(self.btnfrm, text="save", command=self.save_spec)
        self.savebtn.pack(side="left")

        fig = Figure(figsize=(2, 1))
        self.canvas = {
            "si": tk.Canvas(self, width=self.pixdim[0], height=self.pixdim[1]),
            "ax-": tk.Canvas(self, width=self.pixdim[0], height=self.pixdim[1]),
            "ax0": tk.Canvas(self, width=self.pixdim[0], height=self.pixdim[1]),
            "ax+": tk.Canvas(self, width=self.pixdim[0], height=self.pixdim[1]),
            "spc": FigureCanvasTkAgg(fig, self),
        }
        self.axes = {"spc": fig.add_subplot(111)}

        self.btnfrm.pack(side="top")

        # if we just use an embeded bind lambda
        # only the last one created will be used for all
        def mklmd(k):
            return lambda e: self.img_click(k, e)

        # special keys only on main axial
        self.canvas["ax0"].bind("<Return>", lambda e: self.inc_roi_selected())
        self.canvas["ax0"].bind("<Up>", lambda e: self.move_roi(y=-1))
        self.canvas["ax0"].bind("<Down>", lambda e: self.move_roi(y=1))
        self.canvas["ax0"].bind("<Left>", lambda e: self.move_roi(x=-1))
        self.canvas["ax0"].bind("<Right>", lambda e: self.move_roi(x=1))

        # add click to all images
        # left to place, right to go to next roi
        for k, c in self.canvas.items():
            if type(c) == tk.Canvas:
                c.bind("<Button-1>", mklmd(k))
                c.bind("<Button-3>", lambda e: self.inc_roi_selected())

        # dropdown for roi selection
        self.roiselect = tk.Listbox(self)
        self.roiselect.bind(
            "<<ListboxSelect>>", lambda e: self.i_curroi.set(e.widget.curselection()[0])
        )
        self.roiselect.pack(side=tk.TOP, expand=1)

        # layout
        self.canvas["ax+"].pack(side=tk.LEFT)
        self.canvas["ax0"].pack(side=tk.LEFT)
        self.canvas["ax-"].pack(side=tk.LEFT)
        self.canvas["si"].pack(side=tk.BOTTOM)
        self.canvas["spc"].get_tk_widget().pack(side=tk.LEFT)

    def read_ni(self):
        """populate self.imgs dictionary with each loaded neuroimage"""
        if self.fnames["t1"]:
            self.t1 = np.rot90(nib.load(self.fnames["t1"]).dataobj)
            self.pixdim = self.t1.shape
            center = self.pixdim[2] / 2
            hlf = self.voxdim[2] / 2
            self.imgs["ax-"] = npimg(self.t1[:, :, int(center - hlf)])
            self.imgs["ax0"] = npimg(self.t1[:, :, int(center)])
            self.imgs["ax+"] = npimg(self.t1[:, :, int(center + hlf)])
        else:
            ax = np.ones((self.pixdim[0], self.pixdim[1])) * 150
            self.imgs["ax-"] = ImageTk.PhotoImage(image=Image.fromarray(ax))
            self.imgs["ax0"] = ImageTk.PhotoImage(image=Image.fromarray(ax))
            self.imgs["ax+"] = ImageTk.PhotoImage(image=Image.fromarray(ax))
        # self.img['sag'] =  ImageTk.PhotoImage(image=Image.fromarray(tmpax))
        # self.img['cor'] =  ImageTk.PhotoImage(image=Image.fromarray(tmpax))

        if self.fnames["si"]:
            self.siarray = SIArray(self.fnames["si"])
            intgrl = self.siarray.integrateSI(0)
            cm = plt.get_cmap("viridis")
            res = cv2.resize(intgrl, self.pixdim[0:2], interpolation=cv2.INTER_NEAREST)
            colored = cm(res / res.max()) * 255
            im = Image.fromarray(colored[:, :, 0:3].astype(np.uint8))

            # things to hold onto
            self.siIntg = res
            self.imgs["si"] = ImageTk.PhotoImage(image=im)
        else:
            ax = np.ones((self.pixdim[0], self.pixdim[1])) * 150
            self.siIntg = None
            self.imgs["si"] = ImageTk.PhotoImage(image=Image.fromarray(ax))

    def draw_images(self):
        """redraw all images"""
        # redraw image
        for k, c in self.canvas.items():
            if k[0:2] in ["ax", "si"]:
                c.delete("ALL")
                c.create_image(
                    self.pixdim[0], self.pixdim[1], anchor="se", image=self.imgs[k]
                )
            elif k in ["spc"]:
                # dont draw images on matplotlibs
                pass
            else:
                print("no draw code for %s" % k)

    def update_roi_label(self):
        "set current roi label to include box position"
        lb = self.roiselect
        i = lb.curselection()

        # update might happen before listbox has any selection
        if not i:
            print(f"WARN: update roi_label but no i!")
            return
        # listbox curselection is (index, None)
        i = i[0]
        title = self.coords[i].label(self.scout.res if self.scout else 216)

        # no way to change label? rm and add back
        # color is cleared with delete, need to restore
        lb.delete(i)
        lb.insert(i, title)
        lb.itemconfig(i, {"bg": self.hexcolors[i]})

    def update_plot(self):
        """update matplotlib objects"""
        # TODO: get pos from clicked loc
        if not self.scout:
            return
        pos = np.array([self.coords[0].xy])
        (spectrums, fnames) = self.siarray.ReconCoordinates3(self.scout, pos)
        self.axes["spc"].clear()
        self.axes["spc"].plot(spectrums[0])
        self.canvas["spc"].draw()

    def add_coords(self):
        """draw boxes for each coordinage"""
        for i, roi in enumerate(self.coords):
            xy = roi.xy

            # white if selected otherwise same as roiselection
            if i == self.i_curroi.get():
                color = "white"
            else:
                color = self.hexcolors[i]

            # boxes on each canvas
            for k, c in self.canvas.items():
                # skip matplotlib objects
                if k in ["spc"]:
                    continue

                c.create_rectangle(
                    xy[0] - self.voxdim[0] / 2,
                    xy[1] - self.voxdim[1] / 2,
                    xy[0] + self.voxdim[0] / 2,
                    xy[1] + self.voxdim[1] / 2,
                    outline=color,
                )

            # circles for independance
            self.canvas["ax0"].create_circle(
                xy[0], xy[1], self.voxdim[0], outline="red"
            )

    def update(self):
        self.draw_images()
        self.add_coords()
        self.update_plot()
        self.roiselect.selection_set(self.i_curroi.get())
        self.update_roi_label()

    def load(self):
        """read imaging files"""
        if not self.fnames["t1"] or not self.fnames["si"]:
            print("WARNING: missing t1 or si. cannot load")
            return
        self.read_ni()
        self.update()
        self.scout = Scout(None, res=self.pixdim[0])  # 216

    def save_spec(self, outdir="out"):
        "write positioned coordinates recon spectrum.xx.yy files"
        pos = np.array([c.sid3(self.scout.res) for c in self.coords])
        (specs, fnames) = self.siarray.ReconCoordinates3(self.scout, pos, outdir)
        print(specs)
        lcmodel.run_lcmodel(fnames)

    def set_coords(self, roixy_list=None):
        if not roixy_list:
            print("WARNING: no rois to show!")
            return

        self.coords = [ROI(label, [x, y]) for (label, x, y) in roixy_list]
        mn = self.roiselect
        mn.delete(0, "end")

        cm = plt.get_cmap("hsv")
        colors = cm(np.arange(len(self.coords)) / len(self.coords)) * 254
        self.hexcolors = [
            "#%02x%02x%02x" % tuple(rgb[0:3]) for rgb in colors.astype(int).tolist()
        ]
        sctres = self.scout.res if self.scout else 216
        for i, roi in enumerate(self.coords):
            mn.insert("end", roi.label(sctres))
            mn.itemconfig(i, {"bg": self.hexcolors[i]})
            # mn.add_command(label=roi, command=setroi(i))

    def img_click(self, canvas_key, event):
        """click moves the roi box"""
        if canvas_key[0:2] in ["ax", "si"]:
            self.coords[self.i_curroi.get()].update(event)
        # event.x,event.y might need to be rearranged for non ax slices

        self.update()
        # set focus so keys (arrows) can be used
        self.canvas[canvas_key].focus_set()


def read_rois(rois_list=[], roi_file=None):
    """read rois. assign default x,y if not given (list)
    and/or read from tab delem file with cols: roi,x,y"""
    xdef, ydef, rois = 0, 0, []
    # rois = [ [roi,xdef+=10,ydef+=10] for roi in pargs.roi_list]
    for roi in rois_list:
        xdef += 5
        ydef += 5
        rois.append([roi, xdef, ydef])

    if roi_file:
        with open(roi_file, "r") as f:
            while l := f.readline():
                roi, x, y = l.split("\t")
                rois.append([roi, float(x), float(y)])
    return rois


def parse_args(args):
    "read arguments for displaying grid"
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-r", "--ref", dest="ref_fname", help="Scout aligned reference image"
    )
    parser.add_argument(
        "-s", "--siarray", dest="si_fname", help="MRSI siarray.1.1 file"
    )
    parser.add_argument(
        "-i",
        "--roi_initial",
        dest="roi_file",
        default=None,
        help="tab delim roi,x,y coordates. x,y are initial guesses",
    )
    parser.add_argument(
        "-l",
        "--rois",
        dest="rois_list",
        nargs="+",
        default=[],
        help="roi labels if no initial guess, alt to --roi_initial",
    )
    parser.add_argument(
        "-g",
        "--gm_mask",
        dest="gm_file",
        default=None,
        help="gray matter mask. same res as reference image",
    )
    parser.add_argument(
        "--sires",
        dest="sires",
        default=24,
        help="resolution of SI (matrix size, symetrical)",
    )

    pargs = parser.parse_args(args)

    rois = read_rois(pargs.rois_list, pargs.roi_file)
    if not rois:
        raise Exception("Must specify --rois or valid/nonempty --roi_initial file")

    pargs.rois = rois

    return pargs


# TODO:
#  * load GM+GM count
#    - toggle gm mask
#  * move to best GM
#  * collision
#  * right click for closest to click, not next num
if __name__ == "__main__":
    import sys

    pargs = parse_args(sys.argv[1:])

    root = tk.Tk()
    app = App(
        master=root, roixy_list=pargs.rois, ref=pargs.ref_fname, si=pargs.si_fname, sires=pargs.sires
    )
    app.mainloop()
