echo 'str_post,cor,cor_mean,overlap'>fit.dat
grep --text 'correlation =' fit.log > fit.logdat2
awk '{print 31-NR","$3 $8 $11}' fit.logdat2 >> fit.dat
rm fit.logdat2

echo 'str_post,cor,cor_mean,overlap'>fit_flipped.dat
grep --text 'correlation =' fit_flipped.log > fit_flipped.logdat2
awk '{print 31-NR","$3 $8 $11}' fit_flipped.logdat2 >> fit_flipped.dat
rm fit_flipped.logdat2

/usr/bin/R --vanilla --slave < plot_fit.r
/usr/bin/R --vanilla --slave < plot_fit_flipped.r