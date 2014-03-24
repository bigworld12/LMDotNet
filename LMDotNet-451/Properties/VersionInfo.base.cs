using System.Reflection;

[assembly: AssemblyVersion("1.7.2.$REVNUM$")]

#if DEBUG
[assembly: AssemblyConfiguration("Dbg; $REVID$; $DATETIME$")]
#else
[assembly: AssemblyConfiguration("Rls; $REVID$; $DATETIME$")]
#endif
