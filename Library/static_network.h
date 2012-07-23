#if !defined(STATIC_static_network_INCLUDED)
#define STATIC_static_network_INCLUDED




class static_network {
	
	
	
	public:
				
		static_network(string);
		static_network(static_network &, int*);
		~static_network();
		
		int connected();
		
		int size() {return dim;};
		double edges() {return tstrength;};
		
		void set_values (bool f) { module_set=f;};
		void set_weighted (bool f) { weighted=f;};
		
		int draw(string, bool);
		int draw(string);
		
		void set_id_value (map <int, int> &);
		
		
			
		
	protected:
	
						
		class  vertex {
				
			public:
				
				vertex(int , int , const deque<int> &, const deque<double> &);
				~vertex();
						
				double kplus(const deque<int> &);
							
				int id_num;
				int value;
				double strength;
				wsarray* links;
				

		};

				
		double kin (const deque<int> &);
		double ktot (const deque<int> &);
		
		int dim;									// number of nodes
		double tstrength;							// number of links (weighted)
		bool weighted;
		bool module_set;
		
		deque <vertex*> vertices;
		
		
			
};





static_network::vertex::vertex(int b, int c, const deque<int> &d, const deque<double> &e) {
	
	id_num=b;
	value=c;
	links=new wsarray(d, e);
	
	strength=0;
	for (int i=0; i<links->size(); i++)
		strength+=links->w(i);
	
}



static_network::vertex::~vertex() {
	
	delete links;
	links=NULL;

}



double static_network::vertex::kplus(const deque<int> &a) {

	// computes the internal degree of the vertex respect with a

	
	double f=0;
	for (int i=0; i<a.size(); i++) {
	
		double w=links->pos_weight_of(a[i]).second;
		if (w>0)
			f+=w;
	}
		
	return f;
	
}




	
static_network::static_network(string file_name) {
	
	
	
	weighted=true;
	module_set=false;
	
	
	int h= file_name.size();
	
	char b[h+1];
	for (int i=0; i<h; i++)
		b[i]=file_name[i];
	b[h]='\0';
	
	
	
	ifstream in(b);

	
	int innum;
	string word;
	
	deque<vector<int> > values;
	
	while(in>>word && word!= "nodes") {}
	
	while(in>>word && word!= "links") {
		
		
		vector<int> first(2);
		in>>innum;
		first[0]=innum;
		in>>innum;
		first[1]=innum;
		
		values.push_back(first);

		
	}
	
	dim=values.size();
	
	tstrength=0;
	
	for (int i=0; i<dim; i++) {
	
		deque<int> row1;
		deque<double> row2;
		
		int one=-2;
		double two=-2;
		
		while (one!=-1) {
			in>> one;
			row1.push_back(one);
		}
		
		while (two!=-1) {
			in>>two;
			row2.push_back(two);
		}
		
		
		
		row1.pop_front();
		row2.pop_front();
		
		row1.pop_back();
		row2.pop_back();
		
		
		
		vertices.push_back(new vertex (values[i][0], values[i][1], row1, row2));
		tstrength+=vertices[i]->strength;
		
		
	}



	tstrength=tstrength/2;
	
//	cout<<"tsrttrta "<<tstrength<<endl;
	
		
}


static_network::static_network(static_network &uno, int* sequence) {
	
	
	
	deque<int> d_;
	deque<double> e_;
	
	
	
	vertex _done_ (0, 0, d_, e_);
	vertices.assign(uno.dim, & _done_);
	
	
	

	weighted=uno.weighted;
	module_set=uno.module_set;

	dim=uno.dim;
	
	
	for (int i=0; i<dim; i++) {
	
		deque<int> row1(uno.vertices[i]->links->size());
		deque<double> row2(uno.vertices[i]->links->size());
		vector < pair <int, double> > sorted_rows(uno.vertices[i]->links->size());
		
		for (int j=0; j<uno.vertices[i]->links->size(); j++)
			sorted_rows[j]= make_pair(   sequence[uno.vertices[i]->links->l(j)]  ,   uno.vertices[i]->links->w(j)  );
		
		
		sort(sorted_rows.begin(), sorted_rows.end());
		
		for (int j=0; j<sorted_rows.size(); j++) {
		
			row1[j]=sorted_rows[j].first;
			row2[j]=sorted_rows[j].second;
		
		}
			
		
		//prints(row1);
		//prints(row2);

		
		
		vertices[sequence[i]]=(new vertex (uno.vertices[i]->id_num, uno.vertices[i]->value, row1, row2));
		tstrength+=uno.vertices[sequence[i]]->strength;
		
		
	}



	tstrength=tstrength/2;
	
		
}



static_network::~static_network() {
	
	for (int i=0; i<vertices.size(); i++) {
		
		delete vertices[i];
		vertices[i]=NULL;
	
	}


}


int static_network::connected () {

	int spannet=0;
	deque <set<int> > partic;
	deque <int> present;
	present.assign(dim,0);

	
	while (spannet!=dim) {
		
		set <int> connected;
		set <int> newcon;
				

		for (int i=0; i<dim; i++)
			if (present[i]==0) {
				connected.insert(i);
				newcon.insert(i);
				present[i]=1;
				break;
			}
			
		
	
		while (newcon.size()!=0) {
			
			set <int> nnewcon=newcon;
			newcon.clear();
			set <int>::iterator it=nnewcon.begin();
			while (it!=nnewcon.end()) {
				
				int near=0; 
				
				while (near!=vertices[*it]->links->size()) {
					present[*it]=1;
					if (connected.insert(vertices[*it]->links->l(near)).second)
						newcon.insert(vertices[*it]->links->l(near));
					near++;
				}
				it++;
			}
		}
		
		partic.push_back(connected);
		spannet+=connected.size();
	}
	
	ofstream con_out("./Library/Files/connected_component.dat");

	con_out<<"connected component = "<<partic.size()<<endl;
	con_out<<"dimensions"<<endl;
	
	
	int max_pcon=0;

	for (int i=0; i<partic.size(); i++) {
		con_out<<partic[i].size()<<"   ";
		if (partic[i].size()>=partic[max_pcon].size())
			max_pcon=i;
	}
	con_out<<endl<<endl;
	
	con_out<<"maximum: "<<max_pcon<<endl<<endl;
	
	con_out<<"nodes"<<endl;
	
	int i=0;
	
	deque <int> ivm;
	set <int>::iterator it=partic[max_pcon].begin();
	
	int cc=0;
	while (it!=partic[max_pcon].end()) {
		ivm.push_back((*it));		// maps old labels into the new ones
		con_out<<cc++<<" "<<vertices[*it]->id_num<<" "<<vertices[*it]->value<<endl;
		it++;
	}
	
	con_out<<"\nlinks"<<endl;
	
	i=0;
	it=partic[max_pcon].begin();
	while (it!=partic[max_pcon].end()) {
		
		
		con_out<<i++<<" ";
		for (int j=0; j<vertices[*it]->links->size(); j++)
			con_out<<lower_bound(ivm.begin(),ivm.end(), vertices[*it]->links->l(j))-ivm.begin()<<" ";
		con_out<<-1<<endl;
		
		con_out<<0<<" ";
		for (int j=0; j<vertices[*it]->links->size(); j++)
			con_out<<vertices[*it]->links->w(j)<<" ";
		con_out<<-1<<endl;

			
			
		it++;
	}

	return 0;
}


double static_network::kin (const deque<int> &set) {
	
	double k=0;
	for (int i=0; i<set.size(); i++)
		k+=vertices[set[i]]->kplus(set);

	return k;

}


double static_network::ktot (const deque<int> &set) {
	
	double k=0;
	for (int i=0; i<set.size(); i++)
		k+=vertices[set[i]]->strength;
	return k;

}


int static_network::draw(string file_name) {

	return draw(file_name, weighted);

}


int static_network::draw(string file_name, bool _weighted_) {
	
	
	
		int h= file_name.size();
		
		char b[h+1];
		for (int i=0; i<h; i++)
			b[i]=file_name[i];
		b[h]='\0';
		
		
		
		ofstream graph_out(b);

		
		if (_weighted_) {
			
			for (int i=0; i<vertices.size(); i++)
				for (int j=0; j<vertices[i]->links->size(); j++)
					graph_out<<vertices[i]->id_num<<" "<<vertices[vertices[i]->links->l(j)]->id_num<<" "<<vertices[i]->links->w(j)<<endl;
		
		}
		
		else {
			
			for (int i=0; i<vertices.size(); i++)
				for (int j=0; j<vertices[i]->links->size(); j++)
					graph_out<<vertices[i]->id_num<<" "<<vertices[(*(vertices[i]->links))[j]]->id_num<<endl;

		
		
		}

	return 0;

}


void static_network::set_id_value (map <int, int> &a) {
	
	for (int i=0; i<dim; i++)
		a.insert(make_pair(vertices[i]->id_num, vertices[i]->value));


}


#endif




