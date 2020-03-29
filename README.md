# Spread of Virus

In this repository, we provide the code for SIR and SEIR models to analyze the spread of a virus. As COVID-19 is all over the word, this code allow you to get powerful information about future scenarios using traditional epidemiological models. We are working in this project with [**@vlandaberry**](https://github.com/vlandaberry).

Recent spread of COVID-19 produces the need to indentify how this propagation will affect a particular country in terms of: how long would it take until the virus desappears, when the maximum number of infected people will take place in time, are policies taken by goverment effective in terms of reducing the velocity of the spread? Is the health system prepared to give response to the demand of their services, and if it is not, how many beds in intensive and intermediate care are needed and when the saturation of the health system will take place?

We also predict the evolution over time of the number of infected and propose a way to adjust the prediction as we have more observations of the numbero of infected in a particular country choosing the parameters of the models that better predict the last observation. 

We provide two codes:
- in SIR.py we provide SIR model (Suceptible, Infected and Recovered model) developed by Kermack y McKendrik
(1927). You can find more information about the theory behind de model in SIR documentation. 
- in SEIR.py we provide SEIR model (Suceptible, Exposed, Infected and Recovered. You can find more information on SEIR documentation. 

We are aware that these models may not been capturing some specificities of the spread of this virus. In particular, in the SIR model, suceptibles get infected by having contact with an infected person and infected are transmiting the virus all the period until they get recovered. This seems not to reflect the COVID-19 if we take into account that in most countries once that one person is get the diagnostic of infected, they are in quarantine and the spread from that person stops at that time. On the other hand, SEIR models consider that there is a latent period of time where the person that get the virus does not infect others, and after this period has passed, they become infected and spread the virus until they are finally recovered. These models adapts to SARS flu, but not quite to COVID-19, because it seems to be contagious without latent period or a small latent period of time, and once the person get the diagnosis, quarantine measures help avoid the spread of the disease. 

In COVID-19 we also have a group of people that does not have sympthoms, so they are not diagnosed as infected but they are indeed spreading the virus.

We are now working in a new model SAIRD, to include these difference. Hope we will make it available soon.




