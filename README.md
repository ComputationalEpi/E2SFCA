# E2SFCA

discount_func.RDS
Time-travel decay function of empirically measured responses to how far people would be willing to travel for the SARS-CoV-2 Vaccine. Responses were limited to [1,90] and only people interesed or considering the vaccine. A monotonic cubic spline was fit to 1 minus the empirical cumulative distribution of these responses. The spline approximates a continuous travel-time delay function (such as the commonly used gaussian function) but captures the drastically different behavioral response to time benchmarks (e.g. 59 minutes vs 61 minutes).

For example:
discount_func(4) returns value 0.978331 suggesting that 97.8% of the population would travel at least 4 minutes to get the vaccine.
