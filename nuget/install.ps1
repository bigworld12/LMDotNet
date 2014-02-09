param($installPath, $toolsPath, $package, $project)

$lmfit = $project.ProjectItems.Item("lmfit32").ProjectItems.Item("lmfit.dll")
$lmfit.Properties.Item("BuildAction").Value = 0 # BuildAction = None
$lmfit.Properties.Item("CopyToOutputDirectory").Value = 1 # CopyToOutputDirectory = Copy always

$lmfit = $project.ProjectItems.Item("lmfit64").ProjectItems.Item("lmfit.dll")
$lmfit.Properties.Item("BuildAction").Value = 0 # BuildAction = None
$lmfit.Properties.Item("CopyToOutputDirectory").Value = 1 # CopyToOutputDirectory = Copy always
