import gurobipy as gp
from gurobipy import GRB

def optimize_test_basic(T, B, Delta):

    print("\n***** Gurobipy solver - Basic*****\n")

    # Declare the model
    
    m = gp.Model("vaccinations")

    # Variables

    first_doses = m.addVars(T, lb=0.0, vtype=GRB.CONTINUOUS, name="First_Doses")
    second_doses = m.addVars(T, lb = 0.0, vtype=GRB.CONTINUOUS, name="Second_Doses")
    stocks = m.addVars(T, lb=0.0, vtype=GRB.CONTINUOUS, name="Stocks")
    
    # Dictionary for the time sloot, needed in the objective function

    dictio = {}
    for i in range(0, T):
        dictio[i] = [i+1]
    
    time_slot, time  = gp.multidict(dictio)

    #       { 0, i >= T - Delta
    # x_i = 
    #       { y_{i+Delta}, i < T - Delta

    m.addConstrs( first_doses[i] == second_doses[i+Delta] for i in range(0, T-Delta))
    m.addConstrs( first_doses[i] == 0 for i in range(T-Delta, T))

    # y_i = { 0, i < Delta

    m.addConstrs( second_doses[i] == 0 for i in range(0, Delta))

    #         { x_0 + s_0 = b_0, i = 0
    # x_i =   { x_i + s_i = b_i + s_{i-1}, 1 <= i < Delta 
    #         { x_i + x_{i-Delta} + s_i = b_i + s_i, i >= Delta

    m.addConstr ( (first_doses[0]  + 0 + stocks [0] == B[0] + 0))
    m.addConstrs( (first_doses[i] + 0 + stocks[i] == B[i] + stocks[i-1] for i in range(1, Delta)))
    m.addConstrs( (first_doses[i] + first_doses[i-Delta] + stocks[i] == B[i] + stocks[i-1] for i in range(Delta, T)))

    # s_T = { s_T = 0

    m.addConstr( (stocks[T-1] == 0) ) 

    # Set the object function

    m.setObjective((second_doses.prod(time)), GRB.MINIMIZE)

    print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        
        print("\n***** Verbose solution printing *****\n")
        
        m.printAttr("X")

        print("\n***** Solution's values list printing *****\n")

        first_doses_values = resultList[:T]
        print("First_doses values: " + str(first_doses_values))

        second_doses_values = resultList[T:2*T]
        print("Second_doses values: " + str(second_doses_values))

        stock_values = resultList[2*T:]
        print("Stocks values: " + str(stock_values))

        object_function_value = m.objVal
        print("Object_function value: " + str(object_function_value))
    
    else:
        print("\n***** No solutions found *****")

def optimize_test_capacity(T, B, Delta, Capacity):

    print("\n***** Gurobipy solver - Capacity*****\n")

    # Declare the model

    m = gp.Model("vaccinations")

    first_doses = m.addVars(T, lb=0.0, vtype=GRB.CONTINUOUS, name="First_Doses")
    second_doses = m.addVars(T, lb = 0.0, vtype=GRB.CONTINUOUS, name="Second_Doses")
    stocks = m.addVars(T, lb=0.0, vtype=GRB.CONTINUOUS, name="Stocks")
    
    # Dictionary for the time sloot, needed in the objective function

    dictio = {}
    for i in range(0, T):
        dictio[i] = [i+1]
    
    time_slot, time  = gp.multidict(dictio)

    # c_i <= x_i + y_i

    m.addConstrs( (first_doses[i] + second_doses[i] <= Capacity[i] for i in range(0, T)))

    #       { 0, i >= T - Delta
    # x_i = 
    #       { y_{i+Delta}, i < T - Delta

    m.addConstrs( (first_doses[i] == second_doses[i+Delta] for i in range(0, T-Delta)))
    m.addConstrs( first_doses[i] == 0 for i in range(T-Delta, T))

    # y_i = { 0, i < Delta

    m.addConstrs( second_doses[i] == 0 for i in range(0, Delta))

    #         { x_0 + s_0 = b_0, i = 0
    # x_i =   { x_i + s_i = b_i + s_{i-1}, 1 <= i < Delta 
    #         { x_i + x_{i-Delta} + s_i = b_i + s_i, i >= Delta

    m.addConstr ( (first_doses[0]  + 0 + stocks [0] == B[0] + 0))
    m.addConstrs( (first_doses[i] + 0 + stocks[i] == B[i] + stocks[i-1] for i in range(1, Delta)))
    m.addConstrs( (first_doses[i] + first_doses[i-Delta] + stocks[i] == B[i] + stocks[i-1] for i in range(Delta, T)))

    # s_T = { s_T = 0

    m.addConstr( (stocks[T-1] == 0) ) 

    m.setObjective((second_doses.prod(time)), GRB.MINIMIZE)

    print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        
        print("\n***** Verbose solution printing *****\n")
        
        m.printAttr("X")

        print("\n***** Solution's values list printing *****\n")

        first_doses_values = resultList[:T]
        print("First_doses values: " + str(first_doses_values))

        second_doses_values = resultList[T:2*T]
        print("Second_doses values: " + str(second_doses_values))

        stock_values = resultList[2*T:]
        print("Stocks values: " + str(stock_values))

        object_function_value = m.objVal
        print("Object_function value: " + str(object_function_value))
    
    else:
        print("\n***** No solutions found *****")

def euristic(T, B, Delta):

    print("\n--------------------\n")
    print("***** Test euristic *****\n")

    second_doses_values = []
    object_function_value = 0

    #       { 0, i < Delta
    # y_i = 
    #       { min (b_{i-Delta}, b_i), i >= Delta

    for i in range(Delta, len(B)):
        second_doses_values.append(min(B[i-Delta], B[i]))

    for j in range(0, len(second_doses_values)):
        object_function_value += second_doses_values[j] * (j+ Delta + 1)
    
    print("Second_doses values: " + str(second_doses_values))

    print("Object function value: " + str(object_function_value))

    print("Unused doses: " + str(sum(B) - sum(second_doses_values)*2))

if __name__ == "__main__":
    optimize_test_basic(5, [40,30,20,10,20], 2)
    optimize_test_capacity(5, [40,30,20,10,20], 2, [10,100,100,100,100])
    euristic(5, [40,30,20,10,20], 2)