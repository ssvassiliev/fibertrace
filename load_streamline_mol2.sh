#!/bin/bash
message='\nload_streamlines.sh -c [diffusion|direction] -s [(float)colorMin] -S [(float)colorMax]\n -f [stream filename]
   The default coloring mehod is by direction\n
   Options -s and -S  override automatic color scaling in VMD\n'

color_choice=direction
color=charge
colorMax=10.0
colorMin=0.0
while getopts hc:s:S:f: opt; do
        case ${opt} in
                h) echo -e $message; exit 1;;
                c) color_choice=${OPTARG};;
                s) colorMin=${OPTARG}; setColorScale=true;;
                S) colorMax=${OPTARG}; setColorScale=true;;
                f) filename=${OPTARG};;

        esac
done

#--------------------------------------------"

xyz_map='proc xyz_map {} {\n  
  set color_start 33\n
  display update off\n    
  for {set z 0} {$z < 10} {incr z} {\n
   for {set y 0} {$y < 10} {incr y} {\n
       for {set x 0} {$x < 10} {incr x} {\n
           set r [expr { $x*0.1 }];  set g [expr { $y*0.1}];  set b [expr {$z*0.1}]\n
           color change rgb [expr 100*$z + 10*$y + $x + $color_start] $r $g $b\n
                 }\n
        }\n
}\n
display update on\n
}\n'

echo -e $xyz_map > script.vmd
echo "display projection orthographic" >> script.vmd
echo "axes location off" >> script.vmd
if [ "$color_choice" == "direction" ]
  then
    echo "xyz_map" >> script.vmd
fi

if [ "$color_choice" == "diffusion" ]
  then
    echo color scale method RWB >> script.vmd
    echo color scale midpoint 1.0 >> script.vmd
    echo color scale min 0.25 >> script.vmd
fi
i=0
 echo "mol new $filename autobonds off" >> script.vmd
 echo "mol delrep 0 top" >> script.vmd
 echo "mol representation Bonds 0.1 10.0" >> script.vmd
 echo "mol addrep top" >> script.vmd
 echo "mol clipplane normal 0 0 $i {0 0 -1}" >> script.vmd
 echo "mol clipplane center 0 0 $i {0 0 65}" >> script.vmd
 echo "mol clipplane status 0 0 $i {1}" >> script.vmd
 echo "mol modcolor 0 $((i++)) $color" >> script.vmd
  if [ $setColorScale ] ; then
    echo "mol scaleminmax top 0 $colorMin $colorMax" >> script.vmd
  fi
echo "mol new roi.pdb" >> script.vmd
echo "mol off $i" >>  script.vmd
echo "mol top $((i/2))" >> script.vmd
echo "scale by 1.55"  >> script.vmd
echo "translate by 0.03 0.1 0" >> script.vmd
#echo "render TachyonInternal vmdscene.tga"  >> script.vmd
#echo "display resetview"  >> script.vmd
vmd -size 700 700 -e script.vmd 
rm script.vmd





