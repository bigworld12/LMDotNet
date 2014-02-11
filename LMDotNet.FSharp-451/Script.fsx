// Learn more about F# at http://fsharp.net. See the 'F# Tutorial' project
// for more guidance on F# programming.

#r @"D:\niemeyer\Documents\Projekte\Pohlig-BoneFit\LMA-for-CSharp\LMDotNet\LMDotNet-451\bin\Release\LMDotNet.dll"
#load "LMA.fs"

open System
open LMDotNet.FSharp

let findMin = { LMA.defaultSettings with verboseOutput = true }
              |> LMA.minimize

fun (p: float[]) (r: float[]) -> r.[0] <- p.[0]
|> findMin [|0.0|] 1

