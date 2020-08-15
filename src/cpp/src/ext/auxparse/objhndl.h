
#ifndef __OBJHNDL_H
#define __OBJHNDL_H

//*************************************************************************

template<class T> class CHandle {
public:

	T* rep;
	int* pcount;

	CHandle() { rep=0; pcount=0;}
	CHandle(T* pp) : rep(pp), pcount(new int) { /* (*pcount)=0; */ (*pcount)=1;}
	CHandle(const CHandle& r) : rep(r.rep), pcount(r.pcount) 
	{ 
		if(pcount != 0) (*pcount)++;
	}

	void destroy()
	{
		if(pcount!=0)
		{
			if(--(*pcount)==0)
			{
				delete rep;
				delete pcount;
				rep=0; pcount=0;
			}
		}
	}

	T* operator->() { return rep;}
	T* obj() { return rep;}

	void bind(const CHandle& r)
	{
		if(rep!=r.rep)
		{
			destroy();
			rep = r.rep;
			pcount = r.pcount;
			(*pcount)++;
		}
	}

	bool isEmpty() { return (rep == 0)? true : false;}

	CHandle& operator=(const CHandle& r) 
	{ 
		bind(r); return *this;
	}

	int operator<(const CHandle& r) { return (rep<r.rep);}
	int operator>(const CHandle& r) { return (rep > r.rep);}

	int operator==(const CHandle& r) { return (rep==r.rep);}
	int operator!=(const CHandle& r) { return (rep != r.rep);}

	~CHandle()
	{
		destroy();
	}
};

//*************************************************************************

template<class T> inline bool operator <(const CHandle<T>& h1, const CHandle<T>& h2)
{
	return (h1.rep < h2.rep); 
}

//*************************************************************************

template<class T> inline bool operator >(const CHandle<T>& h1, const CHandle<T>& h2)
{
	return (h1.rep > h2.rep); 
}

//*************************************************************************

template<class T> inline bool operator ==(const CHandle<T>& h1, const CHandle<T>& h2)
{
	return (h1.rep == h2.rep); 
}

//*************************************************************************

template<class T> inline bool operator !=(const CHandle<T>& h1, const CHandle<T>& h2)
{
	return (h1.rep != h2.rep); 
}

//*************************************************************************

#endif
