#pragma once

struct comparator
{
	unsigned long int first;
	unsigned long int second;
	comparator *next;
};

class b_sort
{
public:
	int count;
	int *tacts;
	comparator *firstComp;
	comparator *currentComp;

	b_sort(int n);
	~b_sort(void);
	void join(unsigned int i0, unsigned int i1, unsigned int shift, unsigned int n0, unsigned int n1);
	void sort(unsigned int i0, unsigned int step, unsigned int n0);
	void push(unsigned int i0, unsigned int i1);
};

