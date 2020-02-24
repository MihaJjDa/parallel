#include <iostream>
#include "b_sort.h"
#include <math.h>
#include <ctime>

using namespace std;

bool checking(int count, b_sort *sortClass)
{
	bool test[25];
	unsigned int mask;
	bool correct = true;

	for (int i = 0; i < (int)pow(2.0, count); ++i) {
		mask = 1;
		for (int j = count - 1; j >=0 ; j--) {
			if (i & mask > 0) {
				test[j] = true;
			} else {
				test[j] = false;
			}
			mask = mask << 1;
		}
		sortClass->currentComp = sortClass->firstComp;
		while (sortClass->currentComp != NULL)
		{
			if (test[sortClass->currentComp->first] != test[sortClass->currentComp->second] &&
				test[sortClass->currentComp->first]) {
				test[sortClass->currentComp->first] = false;
				test[sortClass->currentComp->second] = true;
			}
			sortClass->currentComp = sortClass->currentComp->next;
		}
		for (int j = 0; j < count -1; j++) {
			if (test[j] != test[j + 1] && test[j]) {
				return false;
			}
		}
	}
	return true;

}

int main(int argc, char* argv[])
{
	unsigned int count, first = 0, step = 1, tacts;
	bool good = true;

	if (argc >= 2) {
		count = atoi(argv[1]);
	} else {
		cout << "No arguments" << endl;
		return -1;
	}

	cout << count << " 0 0" << endl; 
	b_sort *sortClass = new b_sort(count);
	sortClass->sort(first, step, count);

	//good = checking(count, sortClass);

	sortClass->currentComp = sortClass->firstComp;
	while (sortClass->currentComp != NULL)
	{
		cout << sortClass->currentComp->first << ' ' << sortClass->currentComp->second << endl;
		sortClass->currentComp = sortClass->currentComp->next;
	}
	tacts = 0;
	for (int i = 0; i < count; ++i) {
		if (sortClass->tacts[i] > tacts)
			tacts = sortClass->tacts[i];
	}
	cout << sortClass->count << "\n" << tacts << endl;
	//cout << good << endl;
    return 0;
}