void load()
{
  gROOT->ProcessLine(".L src/HiRA_BadMap.cpp+");
  gROOT->ProcessLine(".L src/HiRA_PosCali.cpp+");
//gROOT->ProcessLine(".L src/HiRA_Particle.cpp+");
  gROOT->ProcessLine(".L src/Exp_RunInfo.cpp+");
  gROOT->ProcessLine(".L src/HiRA_CheckHitPattern.cpp+");
  gROOT->ProcessLine(".L src/HiRA_GeoEff.cpp+");
  gROOT->ProcessLine(".L src/HiRA_ReactionLost.cpp+");
}
