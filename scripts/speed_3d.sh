mkdir 2d &>/dev/null
stride=`cat config.conf|grep ^stride|awk '{print $2}'`
dt=`cat config.conf|grep ^dt|awk '{print $2}'`
pdb=`cat config.conf|grep ^coordinates_xyz|awk '{print $2}'`
for dcd in dcd/*.dcd
    do
        ~/Desktop/mt/gitmt/scripts/disas_speed/3d22d ${pdb} ${dcd} ${dcd}_ang 2d/${dcd##*/}
    done
rm raw_speed &> /dev/null
for dcd in 2d/*.dcd
    do
        ~/Desktop/mt/gitmt/scripts/disas_speed/disc ${pdb} ${dcd} timeline #| tail -n 1 | awk '{print $1" "$2}' >> raw_speed
    done
for dat in 2d/*.dat
	do 
		mv ${dat} .
	done
rm -r 2d

# n=`wc -l raw_speed|awk '{print $1}'`
# string="scale=20;($(cat raw_speed|sed -e "s/ /\//g"|paste -sd+))/(${n}*10^-12*${stride}*${dt})"
# bc <<< ${string}