# ./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=5 --Xt=0 | tee scale_high_TB5.dat

start="91.1876"
stop=100000
steps=20

TB=5
Xt_values="0 2.44949"
AS_values="0.1184 0.2368 0.4736"
MT_values="173.34 91.1876"
MTmethod_values="0 1"
M3factor_values="1 0.5"

for Xt in $Xt_values ; do
    for AS in $AS_values ; do
        for MT in $MT_values ; do
            for MTmethod in $MTmethod_values ; do
                for M3factor in $M3factor_values ; do
                    output_file="Mh_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"

                    echo "generating $output_file ..."

                    ./scan.sh --parameter=MS --start="$start" \
                              --stop="$stop" --steps=$steps --step-size=log \
                              --TB="$TB" \
                              --Xt="$Xt" \
                              --AS="$AS" \
                              --MT="$MT" \
                              --MTmethod="$MTmethod" \
                              --M3factor="$M3factor" \
                              > "$output_file"

                    echo "plotting ..."

                    gnuplot -e "filename=\"$output_file\"" plot-Mh-MS.gnuplot

                    output_file_yt="Yt_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"

                    echo "generating $output_file ..."

                    ./scan.sh --parameter=MS --start="$start" \
                              --stop="$stop" --steps=$steps --step-size=log \
                              --TB="$TB" \
                              --Xt="$Xt" \
                              --AS="$AS" \
                              --MT="$MT" \
                              --MTmethod="$MTmethod" \
                              --M3factor="$M3factor" \
                              --output=Yu-3:3 \
                              > "$output_file_yt"

                    output_file_vu="vu_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"

                    echo "generating $output_file ..."

                    ./scan.sh --parameter=MS --start="$start" \
                              --stop="$stop" --steps=$steps --step-size=log \
                              --TB="$TB" \
                              --Xt="$Xt" \
                              --AS="$AS" \
                              --MT="$MT" \
                              --MTmethod="$MTmethod" \
                              --M3factor="$M3factor" \
                              --output=HMIX-103 \
                              > "$output_file_vu"

                    output_file_v="v_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"

                    echo "generating $output_file ..."

                    ./scan.sh --parameter=MS --start="$start" \
                              --stop="$stop" --steps=$steps --step-size=log \
                              --TB="$TB" \
                              --Xt="$Xt" \
                              --AS="$AS" \
                              --MT="$MT" \
                              --MTmethod="$MTmethod" \
                              --M3factor="$M3factor" \
                              --output=HMIX-3 \
                              > "$output_file_v"

                    echo "combining files for Yt, vu and v ..."
                    output_file_combined="yt_vu_v_MS_TB-${TB}_Xt-${Xt}_AS-${AS}_MT-${MT}_MTmethod-${MTmethod}_M3factor-${M3factor}.dat"
                    paste "$output_file_yt" "$output_file_vu" "$output_file_v" > "$output_file_combined"

                    echo "plotting ..."

                    gnuplot -e "filename=\"$output_file_combined\"" plot-mt-MS.gnuplot
                done
            done
        done
    done
done
