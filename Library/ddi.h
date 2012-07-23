

template <typename type>
class ddi {
	
	public :
	
	ddi(const double &a, const double &b) : first(a), second(b) {};
	ddi(const double &a, const double &b, const type &sq) : first(a), second(b), third(sq) {};

	~ddi(){};
	
	
	friend ostream& operator << (ostream& output, const ddi &p) {
		
		output<<endl<<p.first<<" "<<p.second<<" "<<p.third;
		return output; 
	};
	
	
	friend bool operator < (const ddi &a, const ddi &b) {
	
	if (a.first<b.first)
		return true;
	else if (a.first==b.first && a.second<b.second)
		return true;
	else
		return false;
		
	};
	
	friend bool operator == (const ddi &a, const ddi &b) {
	
	return (a.first==b.first && a.second==b.second);
	
	};
	
	double first;
	double second;
	type third;
	
		
};






