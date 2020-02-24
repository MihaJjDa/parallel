#include "b_sort.h"
#include <iostream>
#include <cmath>

using namespace std;

b_sort::b_sort(int n)
{
	count = 0;
	firstComp = NULL;
	currentComp = NULL;
	tacts = new int [n];
	for (int i = 0; i < n; ++i) 
		tacts[i] = 0;
}

b_sort::~b_sort(void)
{
	while (firstComp != NULL) 
	{
		comparator *temp = firstComp->next;
		delete firstComp;
		firstComp = temp;
	}
	delete [] tacts;
}

void b_sort::join(unsigned int i0, unsigned int i1, unsigned int shift, unsigned int n0, unsigned int n1) 
{
	unsigned int i;

	if (n0 * n1 < 1)
		return;

	if (n0 == 1 && n1 == 1) {
		push(i0, i1);
		tacts[i0] = tacts[i1] = max(tacts[i0], tacts[i1]) + 1;
		return;
	}

	join(i0, i1, shift * 2, n0 - (n0) / 2, n1 - (n1) / 2);
	join(i0 + shift, i1 + shift, shift * 2, (n0) / 2, (n1) / 2);
	for (i = 1; i < (n0 - 1); i+=2) {
		push(i0+shift*i, i0+shift*(i+1));
		tacts[i0+shift*i] = tacts[i0 + shift*(i+1)] = max(tacts[i0+shift*i], tacts[i0+shift*(i+1)]) + 1;
	}

	if (n0%2==0) {
		push(i0+shift*(n0-1), i1);
		tacts[i0+shift*(n0-1)] = tacts[i1] = max(tacts[i0+shift*(n0-1)], tacts[i1]) + 1;

		i = 1;
	} else {
		i = 0;
	}

	for (;i<n1-1;i+=2){
		push(i1+shift*i, i1+shift*(i+1));
		tacts[i1 + shift*i] = tacts[i1 + shift*(i+1)] = max(tacts[i1+shift*i], tacts[i1+shift*(i+1)]) + 1;
	}

	return;
}


void b_sort::sort(unsigned int i, unsigned int step, unsigned int n) 
{
	if (n < 2)
		return;

	sort(i, step, n / 2);
	sort(i + step * (n / 2), step, n - n / 2);
	join(i, i + step * (n / 2), step, n / 2, n - n / 2);
	return;
}


void b_sort::push(unsigned int i0, unsigned int i1) 
{
	if (firstComp == NULL) {
		firstComp = new comparator;
		firstComp->first = i0;
		firstComp->second = i1;
		firstComp->next = NULL;
		currentComp = firstComp;
	}
	else {
		currentComp->next = new comparator;
		currentComp = currentComp->next;
		currentComp->first = i0;
		currentComp->second = i1;
		currentComp->next = NULL;
	}
	count++;
}

