#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "Matrix.H"
#include "RateMatrix.H"

RateMatrix::RateMatrix()
{
}

RateMatrix::~RateMatrix()
{
}

int 
RateMatrix::initMatrices(double a, double b)
{
	Q=new Matrix(2,2);

	//We assume that 0 is H+, 1 is H-, 2 is no change. 
	gsl_matrix* q = gsl_matrix_alloc(2, 2);
	
	Q->setValue(-1*(a),0,0);
	Q->setValue(a,0,1);

	Q->setValue(b,1,0);
	Q->setValue(-1*b,1,1);

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			gsl_matrix_set(q,i,j,Q->getValue(i,j));
		}
	}
	gsl_vector_complex *eval = gsl_vector_complex_alloc (2);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2, 2);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (2);
       
        gsl_eigen_nonsymmv (q,eval, evec, w);
     
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);	

	V=new Matrix(2,2);
	L=new Matrix(1,2);
	for(int i=0;i<2;i++)
	{
		double eval_i = GSL_REAL(gsl_vector_complex_get (eval, i));
		L->setValue(eval_i,0,i);
		for (int j=0;j<2;j++)
		{
			double eval_ij=GSL_REAL(gsl_matrix_complex_get(evec,i,j));
			V->setValue(eval_ij,i,j);
		}
	}
        gsl_eigen_nonsymmv_free (w);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_matrix_free(q);
	
	cout << "Q Inits" << endl;
	Q->showMatrix();
	cout << "V Inits" << endl;
	V->showMatrix();
	Vinv=V->invMatrix();
	cout << "Vinv Inits" << endl;
	Vinv->showMatrix();
	cout << "L Inits" << endl;
	L->showMatrix();
	
	a_param=a;
	b_param=b;
	return 0;
}

int 
RateMatrix::setPis(vector<double>& vals)
{
	for(int s=0;s<vals.size();s++)
	{
		pis.push_back(vals[s]);
	}
	return 0;
}

int
RateMatrix::setPi(int k,double v)
{
	pis[k]=v;
}

Matrix* 
RateMatrix::getQ()
{
	return Q;
}

Matrix* 
RateMatrix::getEigenVectors()
{
	return V;
}

Matrix*
RateMatrix::getInvEigenVectors()
{
	return Vinv;
}

Matrix* 
RateMatrix::getEigenValues()
{
	return L;
}

int
RateMatrix::estimateQ(map<string,RateMatrix::SuffStat*>& stats)
{
	Matrix allSum(Q->getRowCnt(),Q->getColCnt());
	for(int k=0;k<Q->getRowCnt();k++)
	{
		for(int l=0;l<Q->getColCnt();l++)
		{
			//Now sum of over the branches in suffstat
			double s=0;
			for(map<string,RateMatrix::SuffStat*>::iterator sIter=stats.begin();sIter!=stats.end();sIter++)
			{
				double aval=0;
				//if(k==l)
				//{
				//	aval=sIter->second->T[k];
				//}
				//else
				//{
					map<int,double>* allvals=NULL;
					RateMatrix::SuffStat* ss=sIter->second;
					if(ss->N.size()==0)
					{	
						continue;
					}
					if(ss->N.find(k)==ss->N.end())
					{
				//		cout <<"Cannot find any entries for " << k  << " for " << sIter->first<< endl;
					}
					else
					{
						allvals=ss->N[k];
						aval=(*allvals)[l];
					}
				//}
				//This sum is the sum of all avals for a given k,l combination.
				s=s+aval;
			}
			allSum.setValue(s,k,l);
		}
	}
	double num_a=allSum.getValue(0,1);
	double den_a=allSum.getValue(0,0);
	double num_b=allSum.getValue(1,0);
	double den_b=allSum.getValue(1,1);
	a_param=num_a/den_a;
	b_param=num_b/den_b;
	cout<<"Rate Est: a=" << a_param << endl;
	cout<<"Rate Est: b=" << b_param << endl;
	Q->setValue(-1*(a_param),0,0);
	Q->setValue(a_param,0,1);

	Q->setValue(b_param,1,0);
	Q->setValue(-1*b_param,1,1);

	cout << "Q Est" << endl;

	Q->showMatrix();
	gsl_matrix* q = gsl_matrix_alloc(2, 2);
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			gsl_matrix_set(q,i,j,Q->getValue(i,j));
		}
	}

	gsl_vector_complex *eval = gsl_vector_complex_alloc (2);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2, 2);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (2);
       
        gsl_eigen_nonsymmv (q,eval, evec, w);
     
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);	

	for(int i=0;i<2;i++)
	{
		double eval_i = GSL_REAL(gsl_vector_complex_get (eval, i));
		L->setValue(eval_i,0,i);
		for (int j=0;j<2;j++)
		{
			double eval_ij=GSL_REAL(gsl_matrix_complex_get(evec,i,j));
			V->setValue(eval_ij,i,j);
		}
	}
        gsl_eigen_nonsymmv_free (w);
	gsl_matrix_free(q);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	cout << "L Est" << endl;
	L->showMatrix();
	cout << "V Est" << endl;
	V->showMatrix();
	Vinv=V->invMatrix();
	cout << "Vinv Est" << endl;
	Vinv->showMatrix();
	return 0;
}

double
RateMatrix::getPi(int index)
{
	return  pis[index];
}

