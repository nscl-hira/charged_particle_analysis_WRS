#include "HiRA_Acceptance.hh"
ClassImp(HiRA_Acceptance);

HiRA_Acceptance::HiRA_Acceptance()
{}

HiRA_Acceptance::~HiRA_Acceptance()
{;}

int GetPid(int zid, int nid)
{
    int mParticleZ[12] = {1,1,1,2,2,2,3,3,3,4,4,4};
    int mParticleN[12] = {0,1,2,1,2,4,3,4,5,3,5,6};
    
    int pid=99;
    for (int i=0;i<12;i++)
    {
        if (zid==mParticleZ[i] && nid==mParticleN[i])
        {
            pid=i;
        }
    }
    if (pid==99){cout<<"no such pid."<<endl;}
    return pid;
}
