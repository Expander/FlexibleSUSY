# ./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=5 --Xt=0 | tee scale_high_TB5.dat

TB=5
Xt_values="0 2.44949"
AS_values="0.1184 0.2368 0.4736"
MT_values="173.34 91.1876"

for Xt in $Xt_values ; do
    for AS in $AS_values ; do
        for MT in $MT_values ; do
            for MTmethod in 0 1 ; do
                for M3factor in 1 0.5 ; do
                    output_file="Mh_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"

                    echo "generating $output_file ..."

                    ./scan.sh --parameter=MS --start="91.1876" \
                              --stop=100000 --steps=20 --step-size=log \
                              --TB="$TB" \
                              --Xt="$Xt" \
                              --AS="$AS" \
                              --MT="$MT" \
                              --MTmethod="$MTmethod" \
                              --M3factor="$M3factor" \
                              > "$output_file"

                    echo "plotting ..."

                    gnuplot -e "filename=\"$output_file\"" plot-Mh-MS.gnuplot
                done
            done
        done
    done
done
