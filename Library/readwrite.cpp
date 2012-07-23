

#if !defined(READWRITE_INCLUDED)
#define READWRITE_INCLUDED
#include "lfm.h"



long seed;

void srand4(void) {
	seed=(long)time(NULL);
}

void srand5(int rank) {
	seed=(long)(rank);
}



double ran2(long *idum) {
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[R2_NTAB];
	double temp;

	if(*idum<=0 || !iy){
		if(-(*idum)<1) *idum=1*(*idum);
		else *idum=-(*idum);
		idum2=(*idum);
		for(j=R2_NTAB+7;j>=0;j--){
			k=(*idum)/R2_IQ1;
			*idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
			if(*idum<0) *idum+=R2_IM1;
			if(j<R2_NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/R2_IQ1;
	*idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
	if(*idum<0) *idum+=R2_IM1;
	k=(idum2)/R2_IQ2;
	idum2=R2_IA2*(idum2-k*R2_IQ2)-k*R2_IR2;
	if (idum2 < 0) idum2 += R2_IM2;
	j=iy/R2_NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum;
	if(iy<1) iy+=R2_IMM1;
	if((temp=R2_AM*iy)>R2_RNMX) return R2_RNMX;
	else return temp;
}



double ran4(void) {
	double r;
	
	r=ran2(&seed);
	return(r);
}


template <typename uno, typename due>
void prints(map <uno, due> &sq,  ostream &out) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		out<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	out<<endl;
	
}



template <typename Seq>
void prints(Seq &sq, ostream &out) {

	typename Seq::iterator it = sq.begin(); 
	while(it != sq.end())
		out<<*(it++)<<"\t";
		

	out<<endl;
	
}

template <typename type_>
void prints(type_ *a, int b) {
	
	for (int i=0; i<b; i++)
		cout<<*(a+i)<<"\t";
	cout<<endl;


}





template<typename T, template<typename> class C>
void printm(C<T>& c, ostream &out) { 
	
	typename C<T>::iterator it = c.begin(); 
	while(it != c.end()) { 
		prints(*it, out);
		it++; 
	} 

	out<<endl;
} 
 
	


template <typename uno, typename due>
void prints(map <uno, due> &sq) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		cout<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	cout<<endl;
	
}

template <typename uno, typename due>
void prints(multimap <uno, due> &sq) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		cout<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	cout<<endl;
	
}



template <typename Seq>
void prints(Seq &sq) {

	typename Seq::iterator it = sq.begin(); 
	while(it != sq.end())
		cout<<*(it++)<<"\t";
		

	cout<<endl;
	
}


template<typename T, template<typename> class C>
void printm(C<T>& c) { 
	
	typename C<T>::iterator it = c.begin(); 
	while(it != c.end()) { 
		prints(*it);
		it++; 
	} 

	cout<<endl;
} 
 



double log_factorial (int num) {
	
	double log_result=0;
	for (int i=1; i<=num; i++)
		log_result+=log(i);
	
	return (log_result);
 
}



double log_combination (int n, int k) {
	
	if (k==0)
		return 0;
	
	if (n<k)
		return -1;
	
	if (n-k<k)
		k=n-k;
	
	double log_c=0;
	for (int i=n-k+1; i<=n; i++)
		log_c+=log(i);
		
	for (int i=1; i<=k; i++)
		log_c-=log(i);
		
	return log_c;
}


double binomial(int n, int x, double p) {		//	returns the binomial distribution, n trials, x successes, p probability

	if (p==0)
		if (x==0)
			return 1;
		else
			return 0;
	
	if (p==1)
		if (x==n)
			return 1;
		else
			return 0;

		
	
	
	double log_b=0;
	log_b+=log_combination(n, x)+x*log(p)+(n-x)*log(1-p);
	return (exp(log_b));
	


}


int use_gnuplot (deque <double> &b) {
	
	ofstream out1("gp.dat");
	out1<<"p \"gp2.dat\" w lp"<<endl;
	ofstream out2("gp2.dat");

	
	for (int i=0; i<b.size(); i++)
		out2<<i<<" "<<b[i]<<endl;
	
	system("gnuplot \"gp.dat\"");
	
	return 0;
		


}

double poisson (int x, double mu) {

	
	return (exp(-mu+x*log(mu)- log_factorial(x)));

}


template <typename type>
int log_histogram (deque <type> &c, ostream & out, int number_of_bins) {		// c is the set od data, min is the lower bound, max is the upper one
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	
	deque <int> hist;
	deque <double> hist2;
	deque <double> bins;
	double step=log(min);
	if (max==min)
		max++;
	
	double bin=(log(max)-log(min))/number_of_bins;		// bin width
	
		

	while (step<=log(max)+2*bin) {
	
		bins.push_back(exp(step));
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	
	for (int i=0; i<c.size(); i++) {
		
		
		int index=bins.size()-1;
		for (int j=0; j<bins.size()-1; j++) if(	(fabs(double(c[i])-bins[j])<1e-7) || (	double(c[i])>bins[j]	&&	double(c[i])<bins[j+1]	)	) { // this could be done in a more efficient way
			
			index=j;
			break;
		
		}
		
		
				
		hist[index]++;
		hist2[index]+=double(c[i]);
		
	}
	
	for (int i=0; i<hist.size()-1; i++) {
		
		double h1= bins[i];
		double h2= bins[i+1];
		double number_of_integer=h2-h1;
		
		//if (fabs(h1 - bins[i])<1e-7)
		//	number_of_integer++;
		
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size()*number_of_integer);
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
		
		

		
		
		
	
	
	}
	
	
	
	return 0;

}



template <typename type>
int histogram (deque <type> &c, ostream & out, int number_of_bins, double b1, double b2) {		// c is the set od data, min is the lower bound, max is the upper one
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	deque <int> hist;
	deque <double> hist2;
	deque <double> bins;
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		bins.push_back(step);
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	
	for (int i=0; i<c.size(); i++) {
		
		
		int index=bins.size()-1;
		for (int j=0; j<bins.size()-1; j++) if(	(fabs(double(c[i])-bins[j])<1e-7) || (	double(c[i])>bins[j]	&&	double(c[i])<bins[j+1]	)	) { // this could be done in a more efficient way
			
			index=j;
			break;
		
		}
		
					
		
				
		hist[index]++;
		hist2[index]+=double(c[i]);
		
	}
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size()*bin);
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}


int irand(int n) {

	return (int(ran4()*(n+1)));
	
}


int shuffle (int *due, const int &dim) {		// it sets due as a random sequence of integers from 0 to dim-1
	
	multimap <int, int> uno;
	for (int i=0; i<dim; i++)
		uno.insert(make_pair(rand(), i));
	

	multimap<int, int>::iterator it;
	
	int h=0;
	for (it=uno.begin(); it!=uno.end(); it++)
		due[h++]=it->second;
	


	return 0;

}



int read_File(bool weighted, bool value_def, const string &file_name, const string &file_name2) {	
	
	// reads the input file, and renames the nodes in order to have a sequence from 0 to N-1
	
	
		
		
	int h= file_name.size();
	
	char b[h+1];
	for (int i=0; i<h; i++)
		b[i]=file_name[i];
	b[h]='\0';	

	
	{	// check 
	
	
		ifstream check_in(b);
		if (!check_in.is_open()) {
			cerr<<"ERROR:\t"<<file_name<<" not found"<<endl;
			return -1;
		}
	
	
	}
	
	
	map <int, int> new_names;		// map from the original names to the new ones

	
	{
		ifstream in(b);
		int a, b;
		double c=1;
		
		
		int label=0;
		while (in>>a) {
			
			in>>b;
			
			if (weighted)
				in>>c;
			
			map<int, int>::iterator itf=new_names.find(a);
			if (itf==new_names.end())
				new_names.insert(make_pair(a, label++));
			
			itf=new_names.find(b);
			if (itf==new_names.end())
				new_names.insert(make_pair(b, label++));
		
		}
	
	}
	
	
	
	
	
	int *values[new_names.size()];
		
	deque< map<int, double> > link_weight;	
	map <int, double> third;
	
	for (int i=0; i<new_names.size(); i++) {
		
		link_weight.push_back(third);

		values[i]=new int[2];
		values[i][0]=0;
		values[i][1]=0;
	
	}
	
	
	{
		ifstream in(b);
		int a, b;
		double c=1;
		
		
		while (in>>a) {		//this loop is to insert links and weights
			
			in>>b;
			
			if (weighted)
				in>>c;
			
			if (a!=b) {
				
				link_weight[new_names[a]].insert(make_pair(new_names[b], c));		
				link_weight[new_names[b]].insert(make_pair(new_names[a], c));
		
			}
		}		
	
	}
	
	
	for (map<int, int>::iterator itp=new_names.begin(); itp!=new_names.end(); itp++)
		values[itp->second][0]=itp->first;
	
	
	if (value_def) {
	
		int h2= file_name2.size();
		
		char b2[h2+1];
		for (int i=0; i<h2; i++)
			b2[i]=file_name2[i];
		
			
		b2[h2]='\0';
		
		
		ifstream in2(b2);
		
		if (in2.is_open()) {
			
			int cc;
			
			while (in2>>cc) {
				
				int ccc;
				in2>> ccc;
				map<int, int>::iterator itf=new_names.find(cc);
				
				if (itf==new_names.end()) {
					cerr<<"ERROR:\t"<<file_name2<<" does not contain the same labels as network.dat"<<endl;
					return -1;
				}
				
				
				values[itf->second][1]=ccc;
		
			}
			
		}
		
		else {
			cerr<<"ERROR:\t"<<file_name2<<" not found"<<endl;
			return -1;
		}
	}

	
	
	ofstream netout("./Library/Files/format_net.dat");
	

	netout<<"nodes"<<endl;
	for (int i=0; i<new_names.size(); i++)
		netout<<i<<" "<<values[i][0]<<" "<<values[i][1]<<endl;
	

	netout<<"\nlinks"<<endl;
	
	
	for (int i=0; i<link_weight.size(); i++) {
		
		
		netout<<i<<" ";
		for (map<int, double>::iterator itm=link_weight[i].begin(); itm!=link_weight[i].end(); itm++)
			netout<<(itm->first)<<" ";
		netout<<-1<<endl;
		
		
		netout<<0<<" ";
		for (map<int, double>::iterator itm=link_weight[i].begin(); itm!=link_weight[i].end(); itm++)
			netout<<(itm->second)<<" ";
		netout<<-1<<endl;
	
	
	
	}
	
	
	for (int i=0; i<new_names.size(); i++) {
		delete[] values[i];
		values[i] = NULL;
	}




	return 0;


}


int read_File(bool weighted, const string &s) {
	
	return read_File(weighted, false, s, s);
}


int read_File(const string &s) {
	
	return read_File(false, false, s, s);
}


int read_File(const string &s, const string &s2) {

	return read_File(false, true, s, s2);
}


int read_File(bool weighted, const string &s, const string &s2) {

	return read_File(weighted, true, s, s2);
}






#endif


