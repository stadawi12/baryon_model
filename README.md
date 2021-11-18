# baryon_model
My phd project: hadrons under extreme conditions. 
This is a model of a baryon in mathematica.

There are 2 main code files that are used for generating results: 
`baryon_code.nb` & `sho_code.nb`
Difference between these two are just the potentials used, 
I decided to separate them into two separate codes because 
it is cleaner that way.
You can run everything in those files and see how they work.

Then there are 3 other files that serve as modules: 
`oscillator_bracket.nb` `mews_and_alphas.nb` `operators.nb`
These are imported before running anything in `baryon_code.nb`
or `sho_code.nb` (first cell of each file imports those modules).
