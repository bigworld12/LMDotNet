using System.Reflection;

[assembly: AssemblyVersion("1.7.1.$REVNUM$")]

#if DEBUG
[assembly: AssemblyConfiguration("Dbg; $REVID$; $DATETIME$")]
#else
[assembly: AssemblyConfiguration("Rls; $REVID$; $DATETIME$")]
#endif
