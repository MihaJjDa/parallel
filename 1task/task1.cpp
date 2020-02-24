#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

using namespace std;

typedef pair<int,int> comp;
typedef vector<int> ints;
typedef vector<comp> comps;
typedef vector<comps> tacts;

void makeTacts(int n, comps c, tacts &t)
{
    ints m(n);
    int fst, snd, num;
    comps tmp;
    for (comps::iterator it = c.begin(); it != c.end(); it++) 
    {
        fst = it->first;
        snd = it->second;
        if (m[fst] > m[snd])
            num = m[fst];
        else
            num = m[snd];
        num++;
        m[fst] = num;
        m[snd] = num;
        if (num > t.size())
            t.push_back(tmp);
        t[num-1].push_back(*it);
    }
    return;
}

void makeSort(int first, int step, int count, comps &c);
void makeMerge(int first1, int first2, int step, int count1, int count2, comps &c);

void makeSort(int first, int step, int count, comps &c)
{
    if (count == 2)
        c.push_back(make_pair(first-1, first+step-1));
    else if (count > 2)
    {
        int count1;

        count1 = count/2;
        makeSort(first, step, count-count1, c);
        makeSort(first+step*(count-count1), step, count1, c);
        makeMerge(first, first+step*(count-count1), step, count-count1, count1, c);
    }
    return;
}

void makeMerge(int first1, int first2, int step, int count1, int count2, comps &c)
{

    if (count1 == 1 && count2 == 1)
   { cout << "QKRQ\n";
        c.push_back(make_pair(first1-1, first2-1));}
    else if (count1 > 0 && count2 > 0)
    {
        int odd1, odd2, i;
 
        odd1 = count1 - count1/2; 
        odd2 = count2 - count2/2;
        cout << "IN ODD\n";
        makeMerge(first1, first2, 2*step, odd1, odd2, c);
        cout << "OUT ODD IN NODD\n";
        makeMerge(first1+step, first2+step, 2*step, count1-odd1, count2-odd2, c); cout << "OUT NODD\n";
    
        for (i = 1; i < count1-1; i += 2) {
            cout << "odd\n";
            c.push_back(make_pair(first1+step*i-1, first1+step*(i+1)-1));
    }
        if (count1 % 2 == 0)
        {
            cout << "QQ\n";
            c.push_back(make_pair(first1+step*(count1-1)-1, first2-1));
            i = 1;
        }
        else
            i = 0;

        for (; i < count2-1; i += 2) 
        { cout << "NODD\n";
            c.push_back(make_pair(first2+step*i-1, first2+step*(i+1)-1));}
    }
    return;
}

void add1(ints &m)
{
    int n = m.size() - 1;
    while (n >= 0 && m[n] != 0)
    {
        m[n] = 0;
        n--; 
    }
    if (n >= 0)
        m[n] = 1;
    return;
}

void sort(ints &m, tacts t)
{
    int tmp, fst, snd;
    for (tacts::iterator p = t.begin(); p != t.end(); p++)
        for (comps::iterator q = p->begin(); q != p->end(); q++)
        {
           fst = q->first;
           snd = q->second;
           if (m[fst] > m[snd])
           {
               tmp = m[fst];
               m[fst] = m[snd];
               m[snd] = tmp;
           }
        }
    return;
}

bool isSorted(ints m)
{
    bool ok = true;
    for (int i = 0; i < m.size()-1 && ok; i++)
        ok = m[i] <= m[i+1];
    return ok;
}

int main(int argc, char** argv)
{
    if (argc > 1)
    {
      int n, i;
      comps c;
      tacts t;
      
      sscanf(argv[1], "%d", &n);
 //     makeSort(1, 1, n, c);

      makeMerge(1, 3, 1, 2, 4, c);
      makeTacts(n, c, t);
      cout << n << " 0 0\n";
      for (i = 0; i < c.size(); i++)
          cout << c[i].first << ' ' << c[i].second << endl;
      cout << c.size() << endl << t.size() << endl;
/*
//вывод компараторов потактно
      cout << endl;
      n = 0;
      for (tacts::iterator p = t.begin(); p != t.end(); p++)
      {
        cout << ++n << endl;
        for (comps::iterator q = p->begin(); q != p->end(); q++)
           cout << q->first << ' ' << q->second << endl;
      }
      cout << endl;*/

/*
//проверка
      for (n = 1; n <= 24; n++)
      {
        ints m, tmp;
        c.clear();
        t.clear();
        makeSort(1, 1, n, c);
        makeTacts(n, c, t);
        for (i = 0; i < n; i++)
            m.push_back(0);
        for (i = 0; i < pow(2, n); i++)
        {
            tmp = m;
            sort(tmp, t);
            if (!isSorted(tmp))
                cout << n << ' ' << i << " FAILED\n";
            add1(m);
        }
      }
*/
    }
    return 0;
}
