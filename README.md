# TVA r-package

##  Introduction

TVA stands for Treatment Variant Aggregation. It designates an algorithm developed by [names] in [paper]. The goal is to aggregate treatment variants to better estimate their effect. The general setting is the following:

You want to study the simultaneous effects of $M$ intervention arms on an outcome $y$. Each intervention arm $j \in [1,M]$ has $R_j$ possible dosages. If $R_j=0$, it means the arm $j$ is not activated. In the end, you have $(R_1+1) \cdot \dots \cdot (R_M+1)$ possible unique policies, where a unique policy is a combination of dosages of all intervention arms.

For example, you want to increase the number of vaccination shots $y$ in a population, and you think of 3 ways to do that : 
- Sending reminders to the population through SMS
- Offering a financial incentive to people that get vaccinated
- Spreading information about vaccination benefits through people in the population network

Moreover, you can either send a small or a large amount of SMS, a small or a large financial incentive, and you can spread the information using random people, trusted people, or trusted people that are also central in the population network. 

In that case, $M=3$ and the first arm is "SMS reminder" with $R_1=2$, second arm is "financial incentive" with $R_2=2$ and third arm is "spreading information" with $R_3=3$. Therefore, you have $K=(2+1) \cdot (2+1) \cdot (3+1) = 36$ possible unique policies.

One possibility is to randomly affect each one the 36 unique policies to a subset of the population, and to perform a regression on the data you get: 

$$y=U\cdot\beta + FE \cdot\gamma + \beta_0 + \epsilon$$

Where $U$ is a dummy matrix indicating the unique policy of each observation, $FE$ is a fixed effects matrix, $\beta_0$ is the intercept and $\beta$ is the effect of each unique policy. One of the issue is that the number of unique policies exponentially increases, meaning you would need a lot of data to obtain an accurate estimation of each unique policy real effect. Another issue is that since the standard deviation of each $\beta_i$ is large, the best $\beta_I$ is very likely to be over-estimated.

**TVA algorithm solves this**. In short, TVA will estimate the marginal effects in order to group unique policies together in a small number of pools that contain unique policies with similar effects, and then perform a regression on the pooled data and adjust the estimation for the winner's curse. 


## Quick guide to use TVA

### Your data

Your data must be structured in a certain (but simple way). Specifically, most TVA functions will have the following arguments :
- `data` : a dataframe containing all your observations with columns for the outcome, the dosage in each arm, the fixed effects and the possible weights.
- `arms` : a vector containing the names of the arms columns in your dataframe.
- `fes` : a vector containing the names of the fixed effect columns in your dataframe (equal to `c()` by default, should be left empty if there are no fixed effects). 
- `w` : the name of the weight column (`NULL` by default, should be left `NULL` if  there are no weights). 

### Other arguments

Other arguments can be asked for some functions of the package:
- `compare_to_zero` : a boolean (`FALSE` by default). If `TRUE`, it means you allow pooling two policies that might have a different set of activated arms. If `FALSE`, it means you never want to pool two policies that have a different set of activated arms. For example, the policy `(1,1)` can be pooled with `(0,1)` if `compare_to_zero` is set to `TRUE`, but it can' t be if it's set to `FALSE`. 
- `estimation_function_name` : a string, can be equal to `pval_MSE` or `pval_OSE`. By default, it is equal to `pval_OSE`. This argument allows you to choose the function that will estimate the support of the marginal effects. In short, `pval_OSE` will be fast and `pval_MSE` will take longer, but will provide more coherent results in some edge cases. 
- `cutoff` : a number between `0` and `1`. At first, we suggest choosing  `pval_OSE` and `cutoff=0.05`.  All the estimation functions need a cutoff value that gives some severity on including or not a marginal effect in the support. If you chose `pval_OSE`, the support is estimated by performing an OLS in the marginal space and by taking marginals that have a p-value below `cutoff`. If you chose `pval_MSE`, the support is estimated by repeating OLS in the marginal space, and removing the marginal with the largest p-value at each step, until the largest p-value is smaller than `cutoff`. 

### How to use TVA

Once the data is structured in the way TVA diggests it, we suggest the user doing the following: 

Run `grid_pval_OSE(data=data,arms=arms,fes=fes,y=y,w=w)` or `plot_pval_OSE(data,arms,fes,y,w,compare_to_zero)` to get an idea of the support size you will obtain for different cutoff values.

Then, run `result = do_TVA(data = data, arms = arms, fes = fes, y = y, cutoff = 5 * 10^(-2), w = w, estimation_function_name = 'pval_OSE')` to test the TVA algorithm for a first time. 

Type :
- `result$marginal_support` to see the marginal_support you obtained
- `result$pools_summary` to get some information about your different pools
- `result$unique_policy` to see all your unique policies and to witch pool they belong to
- `result$data` to get you original dataframe but with some new columns giving pooling information
- `result$fes_support` to see the fixed effects  that were estimated to be in the support
- `result$pooled_ols` to get the final OLS on the pooled data
- `result$winners_effect` to get the estimate of the best pool effect and its 95% confidence interval

You can then increase (resp. decrease) the `cutoff` if you obtained too little (resp. too many) pools.



## TVA steps

TVA algorithm can be decomposed in several steps:
1. Transform the data, from unique policies to marginal policies
2. Estimate the non-zero marginal effects
3. Deduce a pooling for our observations
4. Perform an OLS on the pooled data
5. Adjust the best pool effects for the winner's curse

Below are the details of each steps.

### 1. First, TVA will start by breaking down all unique policies into marginal policies. 

The unique policies exist in a unique policy space of dimension $K=[0,R_1] \times \dots \times [0,R_M]$. One can map this unique policy space to a marginal policy space of the same dimension.  Let us specify what we mean by marginal policy space.

A unique policy can be written as $r=(r_1,\dots ,r_M)$. To write this unique policy in the marginal space, we consider all the marginal policies $m = (m_1,\dots ,m_M)$ that $r$ *dominates* ( $r \geq m$, meaning $r_i \geq m_i, \forall i\in[1,M]$). In that case, we say that $m$ *influences* $r$. The decomposition of $r$ in the marginal space is just the sum of all $m$ that *influence* $r$. 

For example, the policy $r(1,2)$ will be decomposed on $m(0,0)+m(0,1)+m(1,0)+m(1,1)+m(1,2)$.

However, in some settings, you might want to only decompose unique policies on marginals that are variants from this unique policy, meaning they have the exact same activated arms. More precisely, you might want to decompose $r$ on the marginals that it *dominates* ($r\geq m$) but which also *resembles* $r$, in the sense of $m_i=0 \Leftrightarrow r_i=0$ 

In that case, in the decomposition of $r(1,2)$, will be $m(1,1) + m(1,2)$.

Wether you want to add this condition or not, one can map a unique policy to a marginals decomposition, and it is easy to see that this mapping is bijective. 

The first step of TVA is therefore to take the input data and create its marginals counterpart, which is represented by the matrix $X$, where there are $K$ columns and $X_{i,j}=1$  if the unique policy $i$ is *influenced* by the marginal $m_j$, and $X_{i,j}=0$ otherwise.

### 2. Secondly, TVA will estimate the marginals effects and determine the support

Once we have represented the data in the marginal space, we will compute the following regression :

$$y=X\cdot\alpha + FE \cdot\gamma + \alpha_0 + \epsilon$$

One can see that with this setup, $\beta_i = \sum\limits_{j \, \textrm{s.t.} \,  m_j \, \textrm{influences} \,  r_i} \alpha_i$. 

We are not interested in the exact values of $\alpha$, but rather the zeros and the non-zeros (more explanations come in step 3). One could actually get the value of $\alpha$ by performing the first OLS on the unique policy space. Therefore, we have two main methods to extract the marginals that have a non-zero $\alpha$:

- The simplest is to compute the OLS regression, and to assume that $\alpha_i=0$ if the corresponding p-value $p_i$ is lower than a cutoff value $p_{cutoff}$. This gives us the list of marginals that have a non-zero effect by comparing the p-values of  the OLS to a choosen cutoff. We call this the *p-value one-step elimination* method.
- A second method, longer but more accurate, is to perform a first OLS, remove the marginal with the highest p-value, perform a second OLS, remove the new marginal with the highest p-value, and repeat until the highest p-value is lower than a cutoff value $p_{cutoff}$. We call this  the *p-value multi-step elimination* method.

At this stage, we now have a list of marginals that have non-zero effects, which we call the support $S_{\alpha}$.


### 3. Then, TVA pools unique policies together based on the support estimation

Once we estimated the support of marginals, we use them to create groups of unique policies. This is the *pooling* stage. 

An intuitive way to pool unique policies is to group the policies that are influenced by the exact same marginals that are in the support. Therefore, if we have two unique policies $r$ and $r'$, then:
$\forall \, m \in S_{\alpha}, m \, \textrm{influences} \, r \Leftrightarrow m \, \textrm{influences} \, r'$

Thus, a pool $P_i$ can be explicitly defined by the subset of $S_{\alpha}$ that influences all the element of $P_i$. Once we fix $S_{\alpha}=(m_{i_1}, \cdots, m_{i_S})$, $P_i$ can be writen as $(x_1, \cdots , x_S)$ where $x_j=1$ if the pool is influenced by $m_{i_j}$, and $0$ otherwise.

At this point, we can compute many characteristics for each pool (number of unique policies by pool, number of observations, minimal unique policy inside the pool etc.)

### 4. TVA can now perform a regression on pooled data

Now that the data is grouped together, we now perform the following regression:

$$y=Z\cdot\eta + FE \cdot\gamma + \eta_0 + \epsilon$$

Where $Z$ indicates the pool of each observation. The confidence interval on each $\eta_i$ will be much smaller than the one of the $\beta_i$ if estimated in the naive way, because we have much more examples by group.

### 5. Finally, TVA adjusts the best effect for the winner's curse

Each pool corresponds to a set of unique policies. At this point, we assume that by implementing any of the unique policies in a pool $i$, the expected effect is $\eta_i$.

The estimated $\hat\eta$ can be used as it is now and if we choose one pool randomly and implement it, we can effectively expect the effect to be $\hat\eta_i$. However, if the user wants to extract the pool $P_I$ with the best effect, then we have a biais coming from the winner's curse. In short, conditional on $P_i$ being the best estimated pool, the probability is high that it happened by chance and that the real effect is below what we estimated (and this probability increases with the number of pool).

For this reason, we implemented the method exposed in [Andrew's paper] and applied it in our use case. Basically, we construct a new estimator of $\eta_I$, using the best estimated effect $\hat\eta_I$, the second best $\hat\eta_J$ and the number of pools. This gives us an unbiased estimator and confidence interval for $\eta_I$ that one has to use if the goal of using TVA is to extract the best pool and implement it.


## Examples

### 1. Construct dummy data

Assume we want to increase the number of vaccination shots $y$ in a population, and we have 3 intervention arms :
- SMS reminders with 2 dosages
- Financial incentive with 2 dosages
- Information compaign with 3 dosages

There are $(2+1)\cdot(2+1)\cdot(3+1)=36$ unique policies.

The ground truth is :
- Any level of SMS has no impact, no matter the levels of the other arms
- Level 1 of financial incentive and no information campaign has a marginal effect of +1
- Level 1 of financial incentive and level 1 of information campaign has a marginal effect of +2.5
- Level 2 of financial incentive and level 1 of information campaign has a marginal effect of +2

Assume we impose the *resemblance* condition, meaning a marginal influences a policy only if they have the exact same activated arms.

Therefore, noting $\beta$ the effect of policy $(a,b,c)$:
- if $1 \leq b$ and $c=0$, then $\beta=1$
- else if $1=b$ and $1 \leq c$, then $\beta=2.5$
- else if $b=2$ and $1 \leq c$, then $\beta=2+2.5=4.5$
- otherwise, $\beta=0$

We simulate $1000$ observations by choosing a random unique policy for each one and $y_i=\beta_i+\epsilon$, where $\epsilon \sim N(0,4)$.
    
Random unique policies:

    n_obs = 1000
    sms = c(0,1,2) %>% sample(.,replace=TRUE, n_obs) 
    incentive = c(0,1,2) %>% sample(.,replace=TRUE, n_obs) 
    information = c(0,1,2,3) %>% sample(.,replace=TRUE, n_obs) 
    df = data.frame(outcome=0, sms=sms, incentive=incentive, information=information )
True outcome:

    df[(df$information == 0) & (df$incentive >=1),'outcome'] = 1
    df[(df$information >= 1) & (df$incentive ==1),'outcome'] = 2.5
    df[(df$information >= 1) & (df$incentive ==2),'outcome'] = 4.5
Add random noise:

    df$outcome = df$outcome + rnorm(n_obs, mean = 0 , sd = 4)

Assume we also have some fixed effects on year and age of the individuals : 

    df$year = c(0,1) %>% sample(.,replace=TRUE, n_obs) 
    df$age  = c(0,1) %>% sample(.,replace=TRUE, n_obs) 
    df[df$year == 1,'outcome'] = df[df$year == 1,'outcome'] + 2
    df[df$age  == 1,'outcome'] = df[df$age  == 1,'outcome'] + 1
    
We can put the data in the unique policy space:

    df$unique_policy_id = as.numeric(factor(paste0(df$sms, df$incentive, df$information)))
    unique_policy_dummies = data.frame(lme4::dummy(df$unique_policy_id))
    df = cbind(df, unique_policy_dummies)

And perform an OLS on those 36 variables and 2 fixed effects:

    unique_policies = names(unique_policy_dummies)
    formula = as.formula(paste0("outcome","~",paste0(c(unique_policies,'year','age'),collapse = "+")))
    OLS = estimatr::lm_robust(formula = formula, data = df)

One will see that it is not straightforward to pool unique policies together with those results, not to mention the fact that the best unique policy effect is largely overestimated.

One can run the TVA package main function :

    arms = c('sms','incentive','information')
    fes = c('year','age')
    y = 'outcome'
    result = do_TVA(data = df, arms = arms, fes = fes, y = y, cutoff = 0.05 )

By default, the support estimation function used is `pval_OSE`, no scaling and no weights.

To check what this produced, one can look at :
- `result$marginal_support` : this gives the list of marginals with non-zero effect and their corresponding `id`
- `result$pools_summary` : this shows, for each pool, the number of unique policies, the number of observation, the list of marginals that influence the pool, the minimal unique policy inside the pool, and some unique policy examples
- `result$unique_policy` : this is the list of all unique policies that were in the data with their corresponding `pool_id`, the list of marginals that influence them and some other information
- `result$data` : this gives back the original dataframe with some new columns. We now have the `pool_id` of each observation, the marginals that influence them, and dummy columns that are indicators of the `pool_id` that were used in the final OLS
- `result$pooled_ols` : this is the final OLS objet, performed on the pooled data
- `result$winners_effect` : this gives the estimate of the best pooled policy effect, corrected for the winner's curse, along with a 95% confidence interval

If not satisfied with the support size, it is of course possible to change the cutoff value (increase it to increase the support size and conversely). To get a rough idea of what one could obtain, it is useful to run the function:

    grid_pval_OSE(data=df,arms=arms,fes=fes,y=y)

It will gives you a grid of cutoffs with their corresponding support size.
It might also be useful to see this visually, with

    plot_pval_OSE(data=df,arms=arms,fes=fes,y=y,scale=FALSE,compare_to_zero=FALSE)
    
As explained previously, the fastest support estimation method is the p-value one-step elimination, but if time is not an issue, it is suggested to use the p-value multi-step elimination (there is no equivalence between the cutoffs that should be used in the p-value one-step elimination and in the multi-step, they have different order of magnitude)

First of all, let's check a few cutoff value and their corresponding support size:

    plot_pval_MSE(data=df,arms=arms,fes=fes,y=y,scale=FALSE,compare_to_zero=FALSE)
    
In my simulation, I see that taking a cutoff of `1e-10` gives a support size of 4. Then, one can call `do_TVA` with those new options:

    result = do_TVA(data = df, arms = arms, fes = fes, y = y, estimation_function_name = 'pval_MSE', cutoff = 1e-10 )

They should give a more accurate estimation of the support than running it with `pval_OSE`.


## Code structure

The main function of the TVA package is `do_TVA`. The arguments it takes are `data`, `arms`, `fes`, `y`, `w`, `cutoff`, `estimation_function_name`, and `compare_to_zero`.

First, `do_TVA` calls `prepare_data`. This function creates an empty marginals matrix by calling `create_empty_marginals_matrix` and pipes it into `fill_marginals_matrix` to put the right values inside it. This matrix is concatenated with the `fes` and `y` columns of the data, an `intercept` column and is then weighted with the `weight_observations` function. 

Then, `do_TVA` estimate the support with either `pval_OSE` or `pval_MSE`. The support is split into a marginals support and a fixed effects support.

The marginals support is used by the `pool_data` function, which assigns each observation to its right pool.

Then, `pools_info` is called to compute some general descriptive information about each pool.

The final OLS is finally computed with `get_pooled_ols`, and the winner's curse algorithm is applied to the best pool effect with `winners_curse`.

In the end, `do_TVA` returns all those objects in a list. 

