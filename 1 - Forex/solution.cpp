#ifndef __PROGTEST__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <algorithm>
#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>
#if defined (__cplusplus) && __cplusplus > 199711L /* C++ 11 */
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <memory>
#include <condition_variable>
#endif /* C++ 11*/
using namespace std;

struct CExchange
{
        CExchange                     ( int               id,
                                        const string    & from,
                                        const string    & to )
         : m_ProblemID ( id ), m_From ( from ), m_To ( to ),
           m_BankIdx ( 0 ), m_Result ( false ), m_Rate ( 0 )
  {}
  void  AddBank                       ( const string    & bank )
  {
    m_Banks . push_back ( bank );
  }
  const int       m_ProblemID;
  const string    m_From;
  const string    m_To;
  vector<string>  m_Banks;

  int             m_BankIdx;
  bool            m_Result;
  vector<string>  m_Currency;
  double          m_Rate;
};
struct CArbitrage
{
                 CArbitrage                    ( int               id,
                                                 const string    & rates )
                  : m_ProblemID ( id ), m_Rates ( rates ), m_Arbitrage ( false ), m_Rate ( 0 )
  {
  }
  const int      m_ProblemID;
  const string   m_Rates;

  bool           m_Arbitrage;
  vector<string> m_Currency;
  double         m_Rate;
};

#endif /* __PROGTEST__ */
#define BUFFER_SIZE 5
#define eps 1e-9

bool equals(double a, double b) {
    return fabs(a - b) < eps;
}


class CGraph
{
  public:
    CGraph(string args) {m_args = args; isArb = false;}
    ~CGraph();
    void LoadGraph();
    void printITable();
    void printPredTable();
    void m_FloydWarshall();
    map<int, string> m_intToStringMap;
    map<string, int> m_stringToIntMap;
    int ** m_predTable;
    double bestRate;
    double arbRate;
    bool isArb;
    int arbIdx;
    double ** m_incidencyTable;

  private:
    string m_args;
    map<string, map<string, double> > m_graph;
    void m_prepareITable(unsigned int tableSize);
};

CGraph::~CGraph()
{
    for(unsigned int i = 0; i < m_stringToIntMap.size(); i++)
        delete [] m_incidencyTable[i];
    delete [] m_incidencyTable;

    for(unsigned int i = 0; i < m_stringToIntMap.size(); i++)
        delete [] m_predTable[i];
    delete [] m_predTable;
}

void CGraph::m_prepareITable(unsigned int tSize)
{
    m_incidencyTable = new double * [tSize];
    for(unsigned int i = 0; i < tSize; i++)
        m_incidencyTable[i] = new double[tSize];

    for(unsigned int i=0; i<tSize; i++)
    {
        for(unsigned int j=0; j<tSize; j++)
            if(i==j) m_incidencyTable[i][j] = 1.0;
            else m_incidencyTable[i][j] = 0.0;
    }

    for(auto it1 = m_graph.begin(); it1 != m_graph.end(); it1++)
    {
        for(auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
        {
            m_incidencyTable[m_stringToIntMap[it1->first]][m_stringToIntMap[it2->first]] = it2->second;
        }
    }


    m_predTable = new int * [tSize];
    for(unsigned int i = 0; i < tSize; i++)
        m_predTable[i] = new int[tSize];


    for(unsigned int i=0; i<tSize; i++)
    {
        for(unsigned int j=0; j<tSize; j++)
        {
            if (/*!equals(m_incidencyTable[i][j],1.0) && */!equals(m_incidencyTable[i][j], 0.0))
            {
                m_predTable[i][j] = j;
            }
            else
            {
                m_predTable[i][j] = -1;
            }
        }
    }

}

void CGraph::m_FloydWarshall()
{
    for(unsigned int k=0; k<m_intToStringMap.size(); k++)
    {
        for(unsigned int i=0; i<m_intToStringMap.size(); i++)
        {
            for(unsigned int j=0; j<m_intToStringMap.size(); j++)
            {
                if (equals(m_incidencyTable[i][k],0.0) || equals(m_incidencyTable[k][j],0.0)) continue;

                if((m_incidencyTable[i][j]+eps) < (m_incidencyTable[i][k] * (m_incidencyTable[k][j])))
				{
					m_incidencyTable[i][j] = m_incidencyTable[i][k] * m_incidencyTable[k][j];
					m_predTable[i][j] = m_predTable[i][k];

					if(i == j && ((m_incidencyTable[i][j] > (1.0001 - eps)) || equals(m_incidencyTable[i][j],(1.0001 + eps))))
					{
                        arbRate = floor(m_incidencyTable[i][j] * 10000)/10000;
                        isArb = true;
                        arbIdx = i;
                        return;
					}
				}
            }
        }
    }
}

void CGraph::printITable()
{
    for(unsigned int i=0; i<m_stringToIntMap.size(); i++)
    {
        for(unsigned int j=0; j<m_stringToIntMap.size(); j++)
            cout << setw(15) << m_incidencyTable[i][j] << "    ";
        cout << endl;
    }
    cout << endl << endl;
}

void CGraph::printPredTable()
{
    for(unsigned int i=0; i<m_stringToIntMap.size(); i++)
    {
        for(unsigned int j=0; j<m_stringToIntMap.size(); j++)
            cout << setw(3) << m_predTable[i][j] << "    ";
        cout << endl;
    }
    cout << endl << endl;
}

struct CProblem
{
    CProblem() {arb = NULL; ex = NULL;}
    CArbitrage * arb;
    CExchange * ex;
};

void CGraph::LoadGraph()
{
    const char * buff = this->m_args.c_str();
    char from[255], to[255];
    double value;
    int s = -1;
    while(true)
    {
        s = -1;
        sscanf(buff, " %[^ -] -> %[^ :] : %lf %*c %n", from, to, &value, &s), buff += s;


        if(m_stringToIntMap.find(from)==m_stringToIntMap.end())
        {
            m_stringToIntMap[from]=m_intToStringMap.size();
            m_intToStringMap[m_intToStringMap.size()] = from;
        }

        if(m_stringToIntMap.find(to)==m_stringToIntMap.end())
        {
            m_stringToIntMap[to]=m_intToStringMap.size();
            m_intToStringMap[m_intToStringMap.size()] = to;
        }

        m_graph[from].insert(make_pair(to,value));


        if(s==-1) break;
    }

    m_prepareITable(m_intToStringMap.size());
}

class CConsultant
{
  public:
    CConsultant  (CArbitrage *(* arbitrageFn)(void), CExchange *(* exchangeFn)(void), void (* completed)( int ));
    static void ExchangeSeq  (CExchange * req);
    static void ArbitrageSeq (CArbitrage * req);
    void Execute (int workers);


    CArbitrage *(* m_arbFn)(void);
    CExchange *(* m_exFn)(void);
    void (* m_completed)( int );
    pthread_t * initThrArr(int workers);
    void destroyThrArr(pthread_t * arr);
    queue<CProblem*> buffer;

    private:
};

CConsultant::CConsultant  (CArbitrage *(* arbitrageFn)(void), CExchange *(* exchangeFn)(void), void (* completed)( int ))
{
    m_arbFn = arbitrageFn;
    m_exFn = exchangeFn;
    m_completed = completed;
}

pthread_t * CConsultant::initThrArr(int workers)
{
    pthread_t * arr = new pthread_t[workers];
    return arr;
}

void CConsultant::destroyThrArr(pthread_t * arr)
{
    delete [] arr;
}

void backtrack(CExchange * req, CGraph * g)
{
    int i=g->m_stringToIntMap[req->m_From];
    int j=g->m_stringToIntMap[req->m_To];

    if(g->m_predTable[i][j]!=-1)
    {
        req->m_Rate = g->m_incidencyTable[i][j];
        req->m_Currency.push_back(req->m_From);
        req->m_Result = true;

        int ind = g->m_stringToIntMap[req->m_From];

        do
        {
            ind = g->m_predTable[ind][j];
            req->m_Currency.push_back(g->m_intToStringMap[ind]);
        } while(ind!=j);
    }
    else
    {
        req->m_Result = false;
    }


}

/*
    int i,j;
    i=g->m_stringToIntMap[req->m_From];
    j=g->m_stringToIntMap[req->m_To];
    if(g->m_predTable[i][j]==-1)
    {
        req->m_Result=false;
        return;
    }
    while(j!=i)
    {
        req->m_Currency.push_back(g->m_intToStringMap[i]);
        i=g->m_predTable[i][j];
    }
    req->m_Currency.push_back(req->m_To);

    for(unsigned int i=0; i<req->m_Currency.size(); i++)
    {
        cout << req->m_Currency[i] << " ";
    }
    cout << endl;*/


void CConsultant::ExchangeSeq(CExchange * req)
{
    vector<CGraph *> graphArr;
    double bestRate = -1;
    int bestBank = -1;
    for(unsigned int i=0; i<req->m_Banks.size(); i++)
    {
        CGraph * g = new CGraph((req->m_Banks[i]));
        graphArr.push_back(g);
    }

    for(unsigned int i=0; i<graphArr.size(); i++)
    {
        graphArr[i]->LoadGraph();

        if(graphArr[i]->m_stringToIntMap.find(req->m_From)==graphArr[i]->m_stringToIntMap.end() || graphArr[i]->m_stringToIntMap.find(req->m_To)==graphArr[i]->m_stringToIntMap.end()) continue;

        graphArr[i]->m_FloydWarshall();
        if(bestRate==0) continue;
        graphArr[i]->bestRate = graphArr[i]->m_incidencyTable[graphArr[i]->m_stringToIntMap[req->m_From]][graphArr[i]->m_stringToIntMap[req->m_To]];
        if(graphArr[i]->bestRate > bestRate)
        {
            bestRate = graphArr[i]->bestRate;
            bestBank = i;
            req->m_Rate = floor(10000*bestRate)/10000;
            req->m_BankIdx = bestBank;
            req->m_Result = true;
        }
    }

    //cout << "Best bank: " << req->m_BankIdx << ",  Best rate " << req->m_Rate << " ";
    if(req->m_Result)
    {
        backtrack(req, graphArr[bestBank]);
    }
    //cout << endl ;
    for(unsigned int i=0; i<req->m_Banks.size(); i++)
    {
        delete graphArr[i];
    }
}

void backtrack(CArbitrage * req, CGraph * g)
{

    if(g->arbIdx != -1)
    {
        int i =g->arbIdx;
        req->m_Currency.push_back(g->m_intToStringMap[i]);

        do
        {
            i = g->m_predTable[i][g->arbIdx];
            req->m_Currency.push_back(g->m_intToStringMap[i]);
        } while(i!=g->arbIdx);
    }
}
/*
    int i,j;
    i=g->arbIdx;
    j=g->arbIdx;
//    g->printPredTable();

    req->m_Currency.push_back(g->m_intToStringMap[i]);
    i=g->m_predTable[i][j];

    while(j!=i)
    {
        req->m_Currency.push_back(g->m_intToStringMap[i]);
        i=g->m_predTable[i][j];
    }

    req->m_Currency.push_back(g->m_intToStringMap[i]);

    for(int i=0; i<req->m_Currency.size(); i++)
        cout << req->m_Currency[i] << " ";
        cout << endl;
*/



void CConsultant::ArbitrageSeq (CArbitrage * req)
{
    CGraph g(req->m_Rates);
    g.LoadGraph();
    g.m_FloydWarshall();
    if(g.isArb)
    {
        req->m_Arbitrage = true;
        req->m_Rate = g.arbRate;
        backtrack(req, &g);
    }
}

queue<CArbitrage*> testQueue;

CArbitrage * testArbProducerFunction()
{
    if(testQueue.size()==0 ) return NULL;
    CArbitrage * ret = testQueue.front();
    testQueue.pop();
    return ret;
}

queue<CExchange*> testExQueue;

CExchange * testExchProducerFunction()
{
    if(testExQueue.size()==0 ) return NULL;
    CExchange * ret = testExQueue.front();
    testExQueue.pop();
    return ret;
}

void completed( int  a)
{
    //cout << "Done: " << a << endl;
}

struct thrArgs
{
    pthread_mutex_t * g_Mtx;
    sem_t * g_Full;
    sem_t * g_Free;
    CConsultant * cons;
};

void * arbProducer(struct thrArgs* arg)

{
    while(true)
    {

        CProblem * p = new CProblem();
        p->ex = NULL;
        p->arb = arg->cons->m_arbFn();
        if(p->arb==NULL)
        {
            delete p;
            return NULL;
        }
        sem_wait ( arg->g_Free );
        pthread_mutex_lock ( arg->g_Mtx );

        arg->cons->buffer.push(p);

        pthread_mutex_unlock ( arg->g_Mtx );
        sem_post ( arg->g_Full );
    }

    return NULL;
}

void exProducer(struct thrArgs* arg)
{
    while(true)
    {
        CProblem * p = new CProblem();
        p->arb = NULL;
        p->ex = arg->cons->m_exFn();
        if(p->ex==NULL)
        {
            delete p;
            return;
        }
        sem_wait ( arg->g_Free );
        pthread_mutex_lock ( arg->g_Mtx );

        arg->cons->buffer.push(p);

        pthread_mutex_unlock ( arg->g_Mtx );
        sem_post ( arg->g_Full );
    }
}

void * consument(struct thrArgs * arg)
{
   while ( 1 )
    {
      CProblem * p;

      sem_wait ( arg->g_Full );
      pthread_mutex_lock ( arg->g_Mtx );
      p = arg->cons->buffer.front();
      if(p->arb==NULL && p->ex==NULL)
      {
        pthread_mutex_unlock ( arg->g_Mtx );
        sem_post ( arg->g_Full );
        return NULL;
      }
      arg->cons->buffer.pop();
      pthread_mutex_unlock ( arg->g_Mtx );
      sem_post ( arg->g_Free );

      if(p->arb==NULL)
      {
        arg->cons->ExchangeSeq(p->ex);
        arg->cons->m_completed(p->ex->m_ProblemID);
      }
      else
      {
        arg->cons->ArbitrageSeq(p->arb);
        arg->cons->m_completed(p->arb->m_ProblemID);
      }
      delete p;

    }

    return NULL;
}


void CConsultant::Execute(int workers)
{
    struct thrArgs arg;
    buffer = queue<CProblem *>();
    pthread_mutex_t g_Mtx;
    sem_t g_Full;
    sem_t g_Free;
    pthread_t * thrArr;

    pthread_attr_t attr;
    pthread_t arbLoader[1];
    thrArr = initThrArr(workers);
    pthread_attr_init ( &attr );
    pthread_attr_setdetachstate ( &attr, PTHREAD_CREATE_JOINABLE );

   pthread_mutex_init ( &g_Mtx, NULL );
   sem_init ( &g_Free, 0, workers*2 );    /* 4 tady byly spatne nastavene limity semaforu */
   sem_init ( &g_Full, 0, 0 );

    arg.cons = this;
    arg.g_Free = &g_Free;
    arg.g_Full = &g_Full;
    arg.g_Mtx = &g_Mtx;



    for ( int i = 0; i < workers; i ++ ) // creating worker threads
    {
        if ( pthread_create ( &thrArr[i], &attr, (void*(*)(void*)) consument, (void*)(struct thrArgs*) &arg ) )
         {
           perror ( "pthread_create" );
           return;
         }
    }

     if ( pthread_create ( &arbLoader[0], &attr, (void*(*)(void*)) arbProducer, (void*)(struct thrArgs*) &arg ) )
     {
       perror ( "pthread_create" );
       return;
     }

//, pthread_mutex_t * g_Mtx, sem_t * g_Full, sem_t * g_Free
    exProducer(&arg);

    pthread_join ( arbLoader[0], NULL );
    CProblem p;
    p.arb = NULL;
    p.ex = NULL;

    sem_wait ( &g_Free );
    pthread_mutex_lock ( &g_Mtx );
    buffer.push(&p);
    pthread_mutex_unlock ( &g_Mtx );
    sem_post ( &g_Full );


    for ( int i = 0; i < workers; i ++ )
        pthread_join ( thrArr[i], NULL );


    pthread_attr_destroy ( &attr );

    pthread_mutex_destroy ( &g_Mtx );
    sem_destroy ( &g_Free );
    sem_destroy ( &g_Full );

    destroyThrArr(thrArr);
}



#ifndef __PROGTEST__
int main ( void )
{

    return 0;
}
#endif /* __PROGTEST__ */
