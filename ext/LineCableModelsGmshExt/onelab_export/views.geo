// Clear only derived field views belonging to this exported project.
// Called by geometry checks and GetDP before each solve, including unchanged runs.
ResultPrefix = StrCat(CurrentDirectory, "results/");
If(PostProcessing.NbViews > 0)
  For view In {PostProcessing.NbViews-1:0:-1}
    If(!StrCmp(StrSub(View[view].FileName,0,StrLen(ResultPrefix)),ResultPrefix))
      Delete View[view];
    EndIf
  EndFor
  Draw;
EndIf
