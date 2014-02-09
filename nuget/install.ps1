param($installPath, $toolsPath, $package, $project)

$lmfit = $project.ProjectItems.Item("native-x86")
$lmfit.Properties.Item("BuildAction").Value = 0 # BuildAction = None
$lmfit.Properties.Item("CopyToOutputDirectory").Value = 1 # CopyToOutputDirectory = Copy always

$lmfit = $project.ProjectItems.Item("native-x64")
$lmfit.Properties.Item("BuildAction").Value = 0 # BuildAction = None
$lmfit.Properties.Item("CopyToOutputDirectory").Value = 1 # CopyToOutputDirectory = Copy always
