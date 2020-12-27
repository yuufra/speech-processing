# echo 'set title "test"' > dat.txt
echo 'set xlabel "freq[Hz]"' > dat.txt
echo 'set ylabel "logarithmic power spectrum"' >> dat.txt
echo 'plot '\''-'\'' with lines linetype 1' >> dat.txt
./a.out speech_sample/A_a.dat 0 512 >> dat.txt
echo 'pause -1' >> dat.txt
gnuplot --persist dat.txt
