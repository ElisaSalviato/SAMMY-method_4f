# INPUT
## Arguments
analysis_name=$1 # Name of the analysis (will be integrated in the plot name)
region_files=$2 # Comma separated list of bed files containing the regions/genes coordinate where to calculated the metaprofile (e.g. gene_list1.bed,gene_list2.bed)
region_labels=$3 # Comma separated list of name that must correspond to the $region_files (e.g. first_list,second_list)
signal_file=$4 # File containg the values used to calculated the metaprofiles (e.g. C002-r1_4f_S2S.bw)
upflk=$5 # Region size in base pair upstream to the central region of the metaprofile (e.g. 2000)
body=$6 # Size in base pair where the regions in which the metaprofile is calculted should "squeeze" (e.g. 3000)
downflk=$7 # Region size in base pair downstream to the central region of the metaprofile (e.g. 2000)
pbin=$8 # Size of the bin in base pair in which $upflk+$body+$downflk shoud be divided (e.g. 10)
cores=$9 # Number of cores will be used in the analysis
colors=$10 # Comma separated list of colors that must correspond to $region_labels (e.g. red,darkred)
ylims=$11 # Comma separated min and max of the Y axis (e.g. 0,10)
plotsizes=$12 # Comma separated length and height of the output plot (e.g. 16,9)
extra=$13 # Extra information that should be added in the title of the plot

### Output variables
out_dir="./Output/"
mat_dir="$out_dir/01-Matrices/"
plot_dir="$out_dir/02-Metaplots/"

### Programs
computeMatrix="/path/to/computeMatrix"
plotProfile="/path/to/plotProfile"

echo "Input data"
echo "Master metafile: $master_metafile"
echo "Analysis row: $nline"
echo "Row: $line"
echo "Ouptut dirs: $out_dir, $mat_dir, $plot_dir"
echo
echo

# PROCESSING
## Make output dir
mkdir -p $out_dir $mat_dir $plot_dir

## Define analysis name
echo "Analysis name: $analysis_name"
echo

## Make the matrix
echo "Matrix info"

### Import the variables needed from the metafile
#### Region file
region_files=$( echo $region_files | tr ',' ' ' )
nregion_files=$( echo $region_files | tr ' ' '\n' | wc -l )
region_labels=$( echo $region_labels | tr ',' ' ' )

echo "Region files: $region_files"
echo "Labels associate to the regions: $region_labels"

#### Signal file
echo "Signal files: $signal_file"

#### Flanking regions and body sizes
echo "Upstream flanking region: $upflk"
echo "Body region: $body"
echo "Downstream flanking region: $downflk"

#### Binsize
echo "Bin: $pbin"

#### Cores used
echo "Number of processors: $cores"

#### Zero skipped
echo "Zeros will be skipped!!!"

### Colors
colors=$( echo $colors | tr ',' ' ' )
echo "Colors: $colors"

### Ylim
ymin=$( echo $ylims | cut -d ',' -f 1 )
ymax=$( echo $ylims | cut -d ',' -f 2 )
echo "Ylims: ${ymin},${ymax}"

### Plot sizes

wid=$( echo $plotsizes | cut -d ',' -f 1 )
hig=$( echo $plotsizes | cut -d ',' -f 2 )
echo "Plot sizes: $wid,$hig"

### Plot title
plot_title="${analysis_name}__${extra}"
echo "Plot title: $plot_title"

### Create the matrix filename
matrix="$mat_dir/${analysis_name}_${pbin}binsize_${body}rbl_${upflk}-${downflk}flkr__${nregion_files}-nregions__$( basename $signal_file | cut -d '.' -f 1 )-signal.mat.gz"
echo "Matrix file: $matrix"
echo

## Compute the matrix
if [[ ! -f "$matrix" ]]
then

    echo "The matrix does not exists, I am creating it!"
    
    $computeMatrix scale-regions -S $signal_file \
    	       -R $region_files \
               --beforeRegionStartLength $upflk \
               --regionBodyLength $body \
               --afterRegionStartLength $downflk \
    	       --numberOfProcessors "max" \
    	       --binSize $pbin \
    	       --skipZeros -o $matrix
    
fi
echo

# OUTPUT
## Make the metaprofile
echo "Making the metaplot"

### Create the matrix filename
meplot="$plot_dir/$( basename $matrix | cut -d '.' -f 1 ).pdf"
echo "Plot file: $meplot"
echo

$plotProfile -m $matrix \
	     -out $meplot \
	     --plotType "lines" \
	     --colors $colors \
	     --plotTitle $plot_title \
	     --regionsLabel $region_labels \
	     --samplesLabel ' ' \
	     --legendLocation "best" \
	     --startLabel "dom start" \
	     --endLabel "dom end" \
	     --yMin $ymin \
	     --yMax $ymax \
	     --plotHeight $hig --plotWidth $wid

