#
# DuffyNGS  ver 1.4.0

VER = 1.4.0

DuffyNGS.tar.gz : 
	rm  -f ./DuffyNGS_*.tar.gz
	#rm  -f src/*.o
	${R_PROGRAM} CMD SHLIB --output=src/DuffyNGS.so src/samtools/*.c  src/*.c 
	${R_PROGRAM} CMD build  --force .
	${R_PROGRAM} CMD INSTALL ./DuffyNGS_${VER}.tar.gz

