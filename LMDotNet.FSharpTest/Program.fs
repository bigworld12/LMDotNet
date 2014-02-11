open System
open LMDotNet.FSharp

[<EntryPoint>]
let main argv = 
    let findMin = { LMA.defaultSettings with verboseOutput = true }
                  |> LMA.minimize

    let res = fun (p: float[]) (r: float[]) -> r.[0] <- p.[1] - p.[0]
                                               r.[1] <- p.[1] - 2.0 * p.[0] + 0.5
              |> findMin [|0.0; 0.0|] 2

    printfn "%s" res.Message
    printfn "%g %g" res.OptimizedParameters.[0] res.OptimizedParameters.[1]
    
    0 // return an integer exit code
