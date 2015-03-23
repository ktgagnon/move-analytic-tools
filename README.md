# move-analytic-tools
A variety of tools used to analyze human or animal movement (e.g., step lengths, turning distributions, etc.).

steps.define and steps.merge can be used to establish step lengths required for testing whether a trajectory approximates a levy walk, brownian motion, a mixture or the two, or something entirely different.

steps.define and steps.merge are based on the "rectangular model" approach reported by:
      Rhee et al. (2011). On the levy-walk nature of human mobility, IEEE/ACM Transactions on Networking, vol. 19,             issue 3.
 
Use the function steps.define first. Then take the list returned from steps.define and feed it to steps.merge. Of course, you can skip the steps.merge function altogether and proceed to fit various distributions to the step lengths.
