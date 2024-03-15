#!/usr/bin/env bash
#
# quick tool tool to generate inputs for grid.py
#
# 20240314WF - init
#
usage(){
    cat <<HEREDOC
    -s --scout     scout anat file w/FOV matching to MRSI acq
    -r --roi       roi atlas (likely MNI)
    -f --fsdir     FreeSurfer subject folder (for skullstripped anat, GM binarize)
    -t --template  template space !!brain-only!! anat matching roi atlas (likely MNI)
    -o --outdir    directory to save outputs to

HEREDOC

}
parseargs(){
  [ $# -eq 0 ] && usage && exit 0

  declare -g scout roi template anat aseg outdir
  # break apart -h and any other single switch option, check '--switch value' pairs exist
  #PARSED=$(getopt\
  #        --options=s:r:f:t:o:h \
  #    --longoptions=scout:,roi:,fsdir:,template:,outdir:,help \
  #    --name "$0" -- "$@")
  #eval set -- "$PARSED"
  while [ $# -gt 0 ]; do
      case "$1" in
          -s|--scout)    scout="$2"; shift 2;;
          -r|--roi)      roi="$2"; shift 2;;
          -f|--fsdir)    fsdir="$2"; shift 2;;
          -t|--template) template="$2"; shift 2;;
          -o|--outdir)   outdir="$2"; shift 2;;
          -h|--help) usage; exit 1;;
          *) echo "unknown option '$1', see $0 -h"; exit 1;;
      esac
  done
  anat="$fsdir/brain.mgz"
  aseg="$fsdir/aseg.mgz"

  roi="${roi:-}"
  for var in scout fsdir anat aseg; do
    [ ! -v $var ] && echo "ERROR: must specify $var" && exit 2
    varval="${!var}"
    [ ! -r "$varval" ] && echo "ERROR: $var '$varval' is not readable/DNE" && exit 2
  done
  [ -n "${roi}" ] && [ ! -r "${roi}" ] &&

  [ -n "$template" -a ! -r "$tempalate" ] &&
      echo "template '$template' does not exist!" && template=""
  if [ -z "$template" ]; then
      echo "WARNING: no/bad ROI template given. try to download MNI 1mm 09c. Will need template flow."
      echo "see      python3 -m pip install templateflow --upgrade --user --break-system-packages"

      # TODO: could use datalad instead?
      # TODO: use resolution of roi file?
      template=$(python3 -c 'import templateflow.api as tf; print(tf.get("MNI152NLin2009cAsym", resolution=[1],desc="brain", suffix="T1w"))')
      #${TEMPLATEFLOW_HOME:-$HOME/.cache/templateflow}/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01-desc-brain_T1w.nii.gz
  fi
  [ -z "$outdir" ] && outdir=$(mktemp -d "${TMP:-/tmp}/MRSIXXXX") && echo "# saving to $outdir"
}
init_output() {
  declare -g scout roi template anat aseg
  mkdir -p "$outdir"
  test -r "$outdir/anat.nii.gz" ||
    mri_convert     "$anat"         "$_"
  test -r "$outdir/gm_mask.nii.gz" ||
    mri_binarize -i "$aseg" --gm -o "$_"
  test -e "$outdir/roi_atlas.nii.gz" ||
    ln -s "$(readlink -f "$roi")"   "$_"
  test -e "$outdir/scout.nii.gz"||
    ln -s "$(readlink -f "$scout")" "$_"
  test -e "$outdir/template.nii.gz" ||
    ln -s "$(readlink -f "$template")" "$_"
}
mrsi_warps(){
  output="${1:?need MRSI align outputdir as first arg}"
  cd "$output"
  # would use 'bet' but dont want to require FSL
  antsBrainExtraction.sh -d 3 -a scout.nii.gz -e brainWithSkullTemplate.nii.gz -o scout_bet.nii.gz


  # bring anat to scout (affine only)
  test -r anat2scout_0GenericAffine.mat ||
    antsRegistrationSyN.sh -d 3 -f scout_bet.nii.gz -m anat.nii.gz -t a -o anat2scout

  antsApplyTransforms  \
    --default-value 0 \
    --input anat.nii.gz \
    --input-image-type 3 \
    --interpolation Linear \
    --output anat_in_scout.nii.gz \
    --reference-image scout.nii.gz \
    --transform anat2scout_0GenericAffine.mat

  ## Bring ROI to scout (2step syn + affine)
  test -r anat2tmpl_1InverseWarp.nii.gz ||
    antsRegistrationSyN.sh -d 3 -f template.nii.gz -m anat.nii.gz -t s -o anat2tmpl
  antsApplyTransforms  \
    --default-value 0 \
    --input roi.nii.gz \
    --input-image-type 3 \
    --interpolation NearestNeighbor \
    --output roi_in_scout.nii.gz \
    --reference-image scout.nii.gz \
    -t anat2scout_0GenericAffine.mat \
    -t "[ anat2tmpl_0GenericAffine.mat, 1 ]" \
    -t anat2tmpl_1InverseWarp.nii.gz \

}

_MRSI_align() {
  declare -g output
  parseargs "$@"
  mrsi_warps "$output"

}

# if not sourced (testing), run as command
if ! [[ "$(caller)" != "0 "* ]]; then
    set -euo pipefail
    trap 'e=$?; [ $e -ne 0 ] && echo "$0 exited in error $e"' EXIT
    _MRSI_align "$@"
    exit $?
fi
