open System
open LMDotNet

[<EntryPoint>]
let main argv = 
    // set up solver
    let solverEngine = { LMA.defaultSettings with verboseOutput = true } 
                       |> LMA.init

    let findMin = LMA.minimize solverEngine
    
    let res1 = fun (p: float[]) (r: float[]) -> r.[0] <- p.[1] - p.[0]
                                                r.[1] <- p.[1] - 2.0 * p.[0] + 0.5
               |> findMin [|0.0; 0.0|] 2

    let res2 = fun (p: float[]) (r: float[]) -> r.[0] <- p.[1] - p.[0]
                                                r.[1] <- p.[1] - 2.0 * p.[0] + 0.5
               |> findMin [|1.0; 1.0|] 2

    printfn "%s" res1.Message
    printfn "%g %g" res1.OptimizedParameters.[0] res1.OptimizedParameters.[1]
    
    0 // return an integer exit code
