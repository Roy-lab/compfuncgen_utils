LFLAG = -lgsl -lgslcblas 
SRC = Framework.C  Gamma.C  Matrix.C  RateEstimator.C  RateMatrix.C  SpeciesDistManager.C

LIBPATH = gsl/lib
INCLPATH1 = gsl/include

CC=g++
CFLAGS = -g

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o estimateRate

