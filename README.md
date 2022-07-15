# IndividualRewrite

rewriting Individual.jl using graph rewriting for rewriting, this repo will be deleted and merged into [Individual.jl](https://slwu89.github.io/Individual.jl/dev/) when all the functionality is recreated, the code has been debugged and tested, and some minimum level of profiling and performance work has been done. 

## To-do

* SIR with birth and death
* SIR on contact graph
* SIR using event scheduling
* For models where events are not conflicting (e.g. simple SIR), check if model can be run without post composing updates to udpate matches.

## Wishlist
    
* each example done with multiple rewriting constructions (SPO,DPO,SqPO)
* either use Catlab to visualize or work out some pretty examples and explanation
* find examples of rewrite rules that require a certain construction (e.g. an SPO rule that canot be a DPO rule, or SqPO rule that cannot be SPO/DPO)
* example where agents are of different types (mosquitoes and humans, or foxes and rabbits)
* examples of migrating/transforming rules
