#!/usr/bin/env python

import tkinter as tk
import numpy as np
from PIL import Image, ImageTk
import nibabel as nib
from siarray import readsi


def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
def _create_box(self, x, y, r, **kwargs):
    return self.create_rectange(x-r, y-r, x+r, y+r, **kwargs)
tk.Canvas.create_circle = _create_circle
tk.Canvas.create_box = _create_box

def npimg(a):
    return ImageTk.PhotoImage(image=Image.fromarray(a))

class App(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        
        # INIT
        self.fnames = {'t1': None, 'gm': None, 'roixy': None}
        self.imgs = {'ax-': None, 'ax0': None, 'ax+': None}
        self.coords = None

        # DIMS
        self.pixdim = (216, 216, 99)
        self.voxdim = (9, 9, 10)
        self.sires = (24,24)
        
        # VIS
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.load = tk.Button(self, text="load", command=self.load)
        
        self.canvas = {
            'ax-': tk.Canvas(self,width=self.pixdim[0],height=self.pixdim[1]),
            'ax0': tk.Canvas(self,width=self.pixdim[0],height=self.pixdim[1]),
            'ax+': tk.Canvas(self,width=self.pixdim[0],height=self.pixdim[1])}

        self.load.pack(side="top")

        # if we just use an embeded bind lambda
        # only the last one created will be used for all
        def mklmd(k):
            return lambda e: self.img_click(k, e) 

        for k,c in self.canvas.items():
            c.bind("<Button-1>", mklmd(k))
            c.pack(side="bottom")
        
        # add more canvases
        self.canvas['si']=tk.Canvas(self,width=self.pixdim[0],height=self.pixdim[1])
        self.canvas['si'].pack(side="left")


    def read_niis(self):
        if self.fnames['t1']:
            self.t1 = np.rot90(nib.load(self.fnames['t1']).dataobj)
            self.pixdim = self.t1.shape
            self.imgs['ax-'] =  npimg(self.t1[:,:,int(self.pixdim[2]/2 - self.voxdim[2]/2)])
            self.imgs['ax0'] =  npimg(self.t1[:,:,int(self.pixdim[2]/2)])
            self.imgs['ax+'] =  npimg(self.t1[:,:,int(self.pixdim[2]/2 + self.voxdim[2]/2)])
        else:
            ax = np.ones((self.pixdim[0],self.pixdim[1]))*150
            self.imgs['ax-'] =  ImageTk.PhotoImage(image=Image.fromarray(ax))
            self.imgs['ax0'] =  ImageTk.PhotoImage(image=Image.fromarray(ax))
            self.imgs['ax+'] =  ImageTk.PhotoImage(image=Image.fromarray(ax))
        #self.img['sag'] =  ImageTk.PhotoImage(image=Image.fromarray(tmpax))
        #self.img['cor'] =  ImageTk.PhotoImage(image=Image.fromarray(tmpax))

    def draw_images(self):
        # redraw image
        for k,c in self.canvas.items():
            if k[0:2] not in ["ax"]:
                continue
            c.delete('ALL')
            c.create_image(self.pixdim[0], self.pixdim[1], anchor="se", image=self.imgs[k])

    def add_coords(self):
        for ar in self.coords:
            # boxes on each canvas
            for c in self.canvas.values():
                c.create_rectangle(ar[0]-self.voxdim[0]/2, ar[1]-self.voxdim[1]/2,
                                   ar[0]+self.voxdim[0]/2, ar[1]+self.voxdim[1]/2,
                                   outline="white")

            # circles for independance 
            self.canvas['ax0'].create_circle(ar[0], ar[1], self.voxdim[0], outline="red")

    def update(self):
        self.draw_images()
        self.add_coords()
    
    def integrateSI(self,s=0, e=None):
        if not e:
           e=self.siarray.shape[0]
        itgr = np.sum(self.siarray[s:e,:],0).reshape(self.sires)

        self.canvas['si'].create_image(self.pixdim[0], self.pixdim[1], anchor="se", image=npimg(itgr))

    def load(self):
        self.fnames['t1'] = "rorig.nii"
        self.fnames['si'] = "siarray.1.1"
        # TODO: load from file (or default to corner)
        self.coords = np.array([[10, 10, 10], [10, 10, 10], [10, 10, 10]])
        self.read_niis()
        # read si array
        self.siarray = readsi(self.fnames['si'])
        self.integrateSI()

        self.update()

    def img_click(self, canvas_key, event):
        if canvas_key[0:2] == 'ax':
            self.coords[0][0] = event.x
            self.coords[0][1] = event.y
        self.update()
        # set focus so keys (arrows) can be used
        self.canvas[canvas_key].focus_set()

# TODO:
#  * load GM+GM count
#  * move to best GM
#  * collision
#  * track box/circle instead of delete/create
# ---
#  * read in siarray
#  * show
#  * spectrum at position

root = tk.Tk()
app = App(master=root)
app.mainloop()
