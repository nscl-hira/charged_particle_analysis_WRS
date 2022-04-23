#include "Exp_RunInfo.hh"
ClassImp(Exp_RunInfo);

Exp_RunInfo::Exp_RunInfo()
{
  SystemNum = 16;
  for(int i=0;i<SystemNum;i++) { RunNum[i] = 0; }
  
  Initial_LabToCM();  
}

Exp_RunInfo::~Exp_RunInfo()
{;}

void Exp_RunInfo::Initial_LabToCM()
{
  
  double Mass_1u = 931.49410242;
  
  const int BeamNum = 2;
  string BeamTag[BeamNum] = {"Ca40","Ca48"};
  double BeamMass[BeamNum] = {39.962590866*Mass_1u,47.95252290*Mass_1u};
  int BeamA[BeamNum] = {40,48};
  
  const int TargetNum = 4; 
  string TargetTag[TargetNum] = {"Ni58","Ni64","Sn112","Sn124"};
  double TargetMass[TargetNum] = {57.935343*Mass_1u,63.9279660*Mass_1u,111.904818*Mass_1u,123.9052739*Mass_1u};
  
  const int BeamEnergyNum = 2; 
  double BeamEnergy[BeamEnergyNum] = {56,140}; //MeV/A
  
  char NameTem[200];
  for(int iBeam=0;iBeam<BeamNum;iBeam++)
  {
    for(int iTarget=0;iTarget<TargetNum;iTarget++)
    {
      for(int iEnergy=0;iEnergy<BeamEnergyNum;iEnergy++)
      {
        int SystemIndex = iBeam*TargetNum*BeamEnergyNum+iTarget*BeamEnergyNum+iEnergy;
        
        sprintf(NameTem,"%s%sE%.0f",BeamTag[iBeam].c_str(),TargetTag[iTarget].c_str(),BeamEnergy[iEnergy]);
        SystemTag_All[SystemIndex] = NameTem;
        BetaZ_LabToCM[SystemIndex] = Cal_BetaZ_LabToCM(BeamMass[iBeam],TargetMass[iTarget],BeamA[iBeam]*BeamEnergy[iEnergy]);
        Rapidity_Beam_Lab[SystemIndex] = Cal_BeamRapidity_Lab(BeamMass[iBeam],BeamA[iBeam]*BeamEnergy[iEnergy]);
        BetaVector_LabToCM[SystemIndex].SetXYZ(0,0,-BetaZ_LabToCM[SystemIndex]);
      }
    }
  }

}

int Exp_RunInfo::GetSystemIndex(string SystemName)
{
  int Index = -1;
  for(int iSys=0;iSys<SystemNum;iSys++)
  {
    if(SystemName==SystemTag_All[iSys]) { Index = iSys; break;}
  }
  if(Index==-1) { cout<<SystemName<<" is not Found in this class, please check again!"<<endl; }
  return Index;
}

double Exp_RunInfo::Cal_BeamRapidity_Lab(double BMass,double BEnergy)
{
  double BeamMomentum = Sqrt(BEnergy*BEnergy+2*BEnergy*BMass);
  TLorentzVector BeamLabTem(0,0,BeamMomentum,BMass+BEnergy);
  double BeamRapidityTem = BeamLabTem.Rapidity();

return BeamRapidityTem;
}

double Exp_RunInfo::Cal_BetaZ_LabToCM(double BMass,double TMass, double BEnergy)
{
  double BeamMomentum = Sqrt(BEnergy*BEnergy+2*BEnergy*BMass);
  double BeamGamma = (BEnergy+BMass)/BMass;
  double Beta = BeamMomentum/(BeamGamma*BMass+TMass);
return Beta;
}

void Exp_RunInfo::AddRunInfo(string SystemTagTem, int StartRun, int EndRun, string BadMap, int ShadowBarTag,string Trigger)
{
  int Index = GetSystemIndex(SystemTagTem);
  for(int i=StartRun;i<=EndRun;i++)
  {
    int NumTem = RunNum[Index];
    for(int ii=0;ii<NumTem;ii++)
    {
      if(RunNO[Index][ii]==i) { cout<<" The "<<SystemTagTem<<" RunNo "<<i<<" already exist load in this class, please check!" <<endl; getchar();}
    }
    
    RunNO[Index][NumTem] = i;
    Run_TriggerCondition[Index][NumTem] = Trigger;
    BadMapVersion[Index][NumTem] = BadMap;
    ShadowBar[Index][NumTem] = ShadowBarTag;
    
    RunNum[Index] = NumTem+1;
  }
}

int Exp_RunInfo::Get_RunNum(string SystemNameTem)
{
  int SysIndex = GetSystemIndex(SystemNameTem);
  return RunNum[SysIndex];
}

int Exp_RunInfo::Get_RunNO(string SystemNameTem,int Index)
{
  int SysIndex = GetSystemIndex(SystemNameTem);
  return RunNO[SysIndex][Index];
}

string Exp_RunInfo::Get_BadMapVersion(string SystemNameTem,int Index)
{
  int SysIndex = GetSystemIndex(SystemNameTem);
  return BadMapVersion[SysIndex][Index];
}

string Exp_RunInfo::Get_Trigger(string SystemNameTem,int Index)
{
  int SysIndex = GetSystemIndex(SystemNameTem);
  return Run_TriggerCondition[SysIndex][Index];
}

void Exp_RunInfo::PrintRunInfo()
{
  for(int iSys=0;iSys<SystemNum;iSys++)
  {
    PrintRunInfo(SystemTag_All[iSys]);
  }
}

void Exp_RunInfo::PrintRunInfo(string SystemNameTem)
{
  int Index = GetSystemIndex(SystemNameTem);
  //Show the dynamic infor:
  cout<<endl<<"****************=>"<<SystemNameTem<<"<=*******************"<<endl;
  cout<<"Index in RunInfo Class: "<<Index<<endl
      <<"Run Number: "<<RunNum[Index]<<endl
      <<"BetaVector_LabToCM : ("<<BetaVector_LabToCM[Index].X()<<","<<BetaVector_LabToCM[Index].Y()<<","<<BetaVector_LabToCM[Index].Z()<<")"<<endl
      <<"Rapidity_Beam_Lab : "<<Rapidity_Beam_Lab[Index]<<endl;
  cout<<"RunInfo : "<<endl;
  
  int PreRunNO = 0;
  string PreBadMap = "Non";
  string PreTrigger = "Non";
  int PreShadowBar = -1;
  bool IsInRunGroup = 0;
  
  for(int i=0;i<RunNum[Index];i++)
  {
    if(BadMapVersion[Index][i]!=PreBadMap || ShadowBar[Index][i]!=PreShadowBar || Run_TriggerCondition[Index][i]!=PreTrigger || RunNO[Index][i]!=PreRunNO+1)
    {
      if(IsInRunGroup==0)
      {
        cout<<"  "<<RunNO[Index][i]<<" - ";
        PreBadMap = BadMapVersion[Index][i];
        PreTrigger = Run_TriggerCondition[Index][i];
        PreRunNO = RunNO[Index][i];
        PreShadowBar = ShadowBar[Index][i];
        IsInRunGroup = 1;
        continue;
      }
      else
      {
        cout<<PreRunNO<<"  "<<PreBadMap<<"  "<<PreShadowBar<<"  "<<PreTrigger<<endl;
        IsInRunGroup = 0;
        
        cout<<"  "<<RunNO[Index][i]<<" - ";
        PreBadMap = BadMapVersion[Index][i];
        PreTrigger = Run_TriggerCondition[Index][i];
        PreRunNO = RunNO[Index][i];
        PreShadowBar = ShadowBar[Index][i];
        IsInRunGroup = 1;
        continue;
      }
    }
    else if(BadMapVersion[Index][i]==PreBadMap && ShadowBar[Index][i]==PreShadowBar && Run_TriggerCondition[Index][i]==PreTrigger && RunNO[Index][i]==PreRunNO+1)
    {
      PreRunNO++;
      IsInRunGroup=1;
      continue;
    }
    
  }
  if(IsInRunGroup==1)
  {
    cout<<PreRunNO<<"  "<<PreBadMap<<"  "<<PreShadowBar<<"  "<<PreTrigger<<endl<<endl;
  }
}

void Exp_RunInfo::Read_RunInfo(string FileWithPath)
{
  //Reset all the run number of each system, this will allow multi-loading of different configure file.
  for(int iRun=0;iRun<20;iRun++) { RunNum[iRun] = 0; }
  
  ifstream infile(FileWithPath);
  cout<<"Reading "<<FileWithPath<<"  ... ";
  
  string System_Tag; int RunNo_Start; int RunNo_End; string BadMap_Tag; int ShadowBar_Tag; string Trigger_Tag;
  
  while(infile)
  {
    infile>>System_Tag;
    //cout<<"   "<<System_Tag<<"  ";
    if(System_Tag=="#") { infile.ignore(99,'\n'); }
    else
    {
      infile>>RunNo_Start>>RunNo_End>>BadMap_Tag>>ShadowBar_Tag>>Trigger_Tag;
      //cout<<RunNo_Start<<"  "<<RunNo_End<<"  "<<BadMap_Tag<<"  "<<ShadowBar_Tag<<"  "<<Trigger_Tag<<endl;
      AddRunInfo(System_Tag, RunNo_Start, RunNo_End, BadMap_Tag, ShadowBar_Tag, Trigger_Tag);
    }
    //the below is only want to make sure the file already point at the end or not.
    infile.ignore(999,'\n');
    char first; infile.get(first);
    infile.putback(first);
  }
  
  cout<<" Finish!"<<endl;
}
