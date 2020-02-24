#include <stdio.h>
#include <vector>
#include <utility>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <float.h>
#include <omp.h>

using namespace std;

typedef pair<int,int> comp;
typedef vector<int> ints;
typedef vector<comp> comps;
typedef vector<comps> tacts;

struct Point 
{
    float coord[2];
    int index;
};

float x(int i, int j)
{
    return (-1) * (i+1) * (j+1);
}

float y(int i, int j)
{
    return (-1) * (i+1) * (j+1);
}


Point* initPoints(int nPer, int rank, int n)
{
    int index = nPer * rank;
    int j;
    Point* points = new Point[nPer];
    
    for (j = 0; j < nPer && index < n; j++)
    {
        points[j].coord[0] = x(index, j);
        points[j].coord[1] = y(index, j);
        points[j].index = index;
        index++;
    }

    for(; j < nPer; j++)
    {
        points[j].coord[0] = FLT_MAX;
        points[j].coord[1] = FLT_MAX;
        points[j].index = -1;
    }
    return points;
}

MPI_Datatype pointType()
{
    MPI_Datatype point;
    MPI_Datatype types[2] = { MPI_FLOAT, MPI_INT };
    int blocks[2] = { 2, 1 };
    MPI_Aint disps[2] = { offsetof(Point, coord), offsetof(Point, index) };
    MPI_Type_create_struct(2, blocks, disps, types, &point);
    MPI_Type_commit(&point);
    return point;
}


void swap(void *pptr1, void *pptr2)
{
    void **ptr1 = (void **) pptr1;
    void **ptr2 = (void **) pptr2;

    void *tmp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = tmp;
}

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
        c.push_back(make_pair(first1-1, first2-1));
    else if (count1 > 0 && count2 > 0)
    {
        int odd1, odd2, i;
 
        odd1 = count1 - count1/2; 
        odd2 = count2 - count2/2; 
        makeMerge(first1, first2, 2*step, odd1, odd2, c);
        makeMerge(first1+step, first2+step, 2*step, count1-odd1, count2-odd2, c);
    
        for (i = 1; i < count1-1; i += 2) 
            c.push_back(make_pair(first1+step*i-1, first1+step*(i+1)-1));
    
        if (count1 % 2 == 0)
        {
            c.push_back(make_pair(first1+step*(count1-1)-1, first2-1));
            i = 1;
        }
        else
            i = 0;

        for (; i < count2-1; i += 2) 
            c.push_back(make_pair(first2+step*i-1, first2+step*(i+1)-1));
    }
    return;
}

void makeHeap(Point *a, Point tmp, int l, int r)
{
    int i = l;
    int j = 2*l+1;  

    while (j <= r)
    {
       if (j < r && a[j].coord[0] < a[j+1].coord[0])
           j++;
       if (tmp.coord[0] < a[j].coord[0])
       {
          a[i] = a[j];
          i = j;
          j = 2 * j + 1;
       }
       else
          break;
    }
    a[i] = tmp;
}

void hsort(Point *a, int n)
{
    int l, r;
    Point tmp;

    l = n / 2 ;
    r = n-1;

    while (l > 0)
    {
        l--;
        tmp = a[l];
        makeHeap(a, tmp, l, r);
    }
    
    while (r >= 0)
    {
        tmp = a[r];
        a[r] = a[0];
        r--;
        makeHeap(a, tmp, 0, r);  
    }
}

void msort(Point *a, Point *b, int n1, int n2)
{
    int i, i1, i2;
    
    Point *tmp = new Point[n1 + n2];
    
    for (i = 0, i1 = 0, i2 = 0; i1 < n1 && i2 < n2;)
    {
        if (a[i1].coord[0] < b[i2].coord[0]) 
        {
            tmp[i] = a[i1];
            i1++;
        }
        else
        {
            tmp[i] = b[i2];
            i2++;
        }
        i++;
    }
    
    for (; i1 < n1; i++, i1++)
        tmp[i] = a[i1];
    for (; i2 < n2; i++, i2++)
        tmp[i] = b[i2];
    
    for (i1 = 0; i1 < n1; i1++)
        a[i1] = tmp[i1];
    for (i2 = 0; i2 < n2; i1++, i2++)
        b[i2] = tmp[i1];
    
    delete[] tmp;
}


void dhsort(Point *a, int n, int threads = 1)
{
    if (n == 2)
    {
        if (a[0].coord[0] > a[1].coord[0])
        {
            Point tmp = a[0];
            a[0] = a[1];
            a[1] = tmp;
        }
    }
    else if (n > 2)
    {
        int count1 = n / 2;
        int count2 = n - count1;
        if (threads == 1)
            if (n > 1e+5)
            {
                dhsort(a, count1);
                dhsort(&(a[count1]), count2);
                msort(a, &(a[count1]), count1, count2);
            }
            else
                hsort(a, n);
        else
        {
            int threads1 = threads/2;
            #pragma omp parallel sections num_threads(threads)
            {
                #pragma omp section
                {
                    dhsort(a, count1, threads1);
                }
                #pragma omp section
                {
                    dhsort(&(a[count1]), count2, threads - threads1);
                }
            }
            msort(a, &(a[count1]), count1, count2);
        }
    }
    return;
}

bool isSorted(Point *points, int n)
{
    bool ok = true;
    for (int i = 1; i < n && ok; i++)
        ok = !(points[i-1].coord[0] > points[i].coord[0]);
    return ok;
}

bool isAllSorted(Point *points, int nPer, int rank, int size, MPI_Datatype MPI_POINT)
{
    Point *left_max = new Point();
    MPI_Sendrecv(points + nPer - 1, 1, MPI_POINT, (rank+1) % size, 1,
                 left_max, 1, MPI_POINT, (rank-1) % size, 1, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int sorted = isSorted(points, nPer) && 
                 (rank == 0 || !(points[0].coord[0] < left_max->coord[0]));
    int numSorted;
    MPI_Allreduce(&sorted, &numSorted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    return numSorted == size;
}


int main(int argc, char **argv)
{
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc > 2 && size > 0)
    {
 
        int n1, n2, threads = 1;
        sscanf(argv[1], "%d", &n1);
        sscanf(argv[2], "%d", &n2);

        if (n1 > 0 && n2 > 0)
        {
            if (argc > 3)
                sscanf(argv[3], "%d", &threads);
            else
                threads = omp_get_max_threads();
            if (threads < 1)
                threads = 1;
            omp_set_num_threads(threads);

            int rank, nPer, i, j, k;
            int n = n1 * n2;
            double startTime, sortTime, endTime;
            double time, dhTime, bTime;
            double maxTime, maxdhTime, maxbTime;
            Point p, p1;
            comps c;
            tacts t;

            MPI_Datatype MPI_POINT = pointType();
            MPI_Status status;

            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            if (n % size > 0)
                nPer = 1 + n / size;
            else
                nPer = n / size;

            Point *points = initPoints(nPer, rank, n);
            Point *tmp = new Point[nPer];
            Point *points1 = new Point[nPer];
/*
            for (int i = 0; i < nPer; i++)
                printf("%d %d %f %f %d\n", 
                       rank, i, points[i].coord[0], points[i].coord[1], points[i].index);*/
                                   
            makeSort(1, 1, size, c);
            makeTacts(size, c, t);
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            startTime = MPI_Wtime();

            dhsort(points, nPer, threads);

            sortTime = MPI_Wtime();

            for (comps::iterator it = c.begin(); it != c.end(); it++) 
            {
                if (rank == it->first) {
                    MPI_Send(points, nPer, MPI_POINT,
                             it->second, 0, MPI_COMM_WORLD);
                    MPI_Recv(points1, nPer, MPI_POINT,
                             it->second, 0, MPI_COMM_WORLD, &status);
                    for (i = 0, j = 0, k = 0; k < nPer; k++) 
                    {
                        p = points[i];
                        p1 = points1[j];
                        if (p.coord[0] < p1.coord[0]) 
                        {
                            tmp[k] = p;
                            i++;
                        } 
                        else 
                        {
                            tmp[k] = p1;
                            j++;
                        }
                    }
                    
                    swap(&points, &tmp);
                }

                if (rank == it->second) {
                    MPI_Recv(points1, nPer, MPI_POINT,
                             it->first, 0, MPI_COMM_WORLD, &status);
                    MPI_Send(points, nPer, MPI_POINT,
                             it->first, 0, MPI_COMM_WORLD);
                    for (i = nPer-1, j = nPer-1, k = nPer-1; k >= 0; k--) 
                    {
                        p = points[i];
                        p1 = points1[j];
                        if (p.coord[0] > p1.coord[0]) 
                        {
                            tmp[k] = p;
                            i--;
                        }  
                        else 
                        {
                            tmp[k] = p1;
                            j--;
                        }
                    }
                    swap(&points, &tmp);
                }
            }

            endTime = MPI_Wtime();
            
            time = endTime - startTime;
            bTime = endTime - sortTime;
            dhTime = sortTime - startTime;

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&dhTime, &maxdhTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&bTime, &maxbTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
/*
            for (int i = 0; i < nPer; i++)
                printf("%d %d %f %f %d\n", 
                       rank, i, points[i].coord[0], points[i].coord[1], points[i].index);*/

            bool sorted = isAllSorted(points, nPer, rank, size, MPI_POINT);
           
            if (rank == 0) 
            {
                cout << "n1:      " << n1        << endl;
                cout << "n2:      " << n2        << endl;
                cout << "n:       " << n         << endl;
                cout << "nPer:    " << nPer      << endl;
                cout << "size:    " << size      << endl;
                cout << "threads: " << threads   << endl;
                cout << "cmps:    " << c.size()  << endl;
                cout << "tacts:   " << t.size()  << endl;
                cout << "dhtime:  " << maxdhTime << endl;
                cout << "btime:   " << maxbTime  << endl;
                cout << "time:    " << maxTime   << endl;
                if (!sorted)
                    cout <<  "not sorted!\n";
            }

            delete [] points;
            delete [] points1;
            delete [] tmp;
        }       
    }

    MPI_Finalize();
    return 0;
}

