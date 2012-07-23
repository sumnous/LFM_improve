


#if !defined(SARRAY_INCLUDED)
#define SARRAY_INCLUDED



#if !defined UNLIKELY 
#define UNLIKELY -2147483647
#endif




// there are two classes here:
// sarray which contains three functions: size(), [], find()
// wsarray which also contains w(), pos_weight_of()



//		about wsarray::
//		the constructor requires two deque (int and double) that must be sorted and have equal size
//		size() gives the size
//		[x] gives the integer in position x (constant time)
//		find(x) gives the position of the integer x, if any; otherwise it gives UNLIKELY (logarithmic time)
//		w(x)	gives the double in position x; l(x) gives the integer in x
//		pos_weight_of(x) gives the position of integer x and its weight





template <typename type>
class sarray {
	
	
	
	public:
				
		sarray(const deque<type> &);
		~sarray();
		
		int size() {return _size_;};
		type operator [] (int i) {			
		
			return _array_[i]; 
		};
		
		type l(int i) {
			
			return _array_[i]; 
		};
		
		int find(type);
		
			
		
	protected:
		
		type* _array_;
		int _size_;
		
};

template <typename type>
sarray<type>::sarray (const deque<type> &a) {
	
	_size_=a.size();
	_array_=new type[_size_];
	
	for (int i=0; i<a.size(); i++)
		_array_[i]=a[i];
	

}



template <typename type>
sarray<type>::~sarray() {
					
	delete[] _array_;
	_array_ = NULL;
				
}

template <typename type>
int sarray<type>::find(type a) {
	
	int one=0;
	int two=_size_-1;
	if (a<_array_[one] || a>_array_[two])
		return UNLIKELY;
	
	if (a==_array_[one])
		return one;
		
	if (a==_array_[two])
		return two;

		
	
	while (two-one>1) {
		
		
		int middle=(two-one)/2 + one;
		
		
		if (a==_array_[middle])
			return middle;

		
		if (a>_array_[middle])
			one=middle;
		else
			two=middle;
	
	}
	
	return UNLIKELY;
				
}



template <typename type>
void prints(sarray<type> &a) {

	for (int i=0; i<a.size(); i++)
		cout<<a[i]<<"\t";
	cout<<endl;

	


} 


template <typename type>
void prints(sarray<type> &a, ostream &out) {

	for (int i=0; i<a.size(); i++)
		out<<a[i]<<"\t";
	out<<endl;

	


}




class wsarray : public sarray<int> {
	
	
	
	public:
				
		wsarray(const deque<int> &, const deque<double > &);
		~wsarray();
		
		double w(int i) {
		
			return _weight_[i];
		
		};
		
		pair <int, double> pos_weight_of (int x);
		
	private:
		
		double* _weight_;
		
};



wsarray::wsarray(const deque<int> &b, const deque<double> &a) : sarray<int>(b) {


	_weight_=new double[sarray<int>::_size_];
	
	for (int i=0; i<a.size(); i++)
		_weight_[i]=a[i];

	

}


wsarray::~wsarray() {
					
	delete[] _weight_;
	_weight_ = NULL;
				
}


pair <int, double>  wsarray::pos_weight_of(int x) {
			
	int i=find(x);
	if (i==UNLIKELY)
		return (make_pair(UNLIKELY, UNLIKELY));
			
	return (make_pair(i, _weight_[i]));
		
		
};



void prints(wsarray &a) {

	for (int i=0; i<a.size(); i++)
		cout<<a[i]<<"\t"<<a.w(i)<<endl;
	cout<<endl;

	


} 


void prints(wsarray &a, ostream &out) {

	for (int i=0; i<a.size(); i++)
		out<<a[i]<<"\t"<<a.w(i)<<endl;
	out<<endl;

	


}




#endif



