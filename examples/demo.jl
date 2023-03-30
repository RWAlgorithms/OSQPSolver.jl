Random.seed!(25)

dic = BSON.load("BLS_prob.bson")
B = dic[:B]
lbs = dic[:lbs]
ubs = dic[:ubs]
y = dic[:observations]