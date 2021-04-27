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
    
    # Dictionary for the time slot, needed in the objective function

    dictio = {}
    for i in range(0, T):
        dictio[i] = i+1
    
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

def optimize_test_capacity_multiple_vaccines(T, B, Delta, Capacity):

    '''
    B = {"Vaccine name": [B_List], "Vaccine name": [B_List] }
    '''

    print("\n***** Gurobipy solver - Capacity Multiple Vaccines*****\n")

    m = gp.Model("vaccinations")

    dicti = {}
    dict_B = {}

    for i in B:
        for j in range(0, T):
            dicti[(i, j)] = [j+1]
            dict_B[(i,j)] = B[i][j]  
    original_B = B
    B = dict_B

    combinations, time_frame = gp.multidict(dicti)

    first_doses = m.addVars(combinations, lb=0.0, vtype=GRB.CONTINUOUS, name="First_Doses")
    second_doses = m.addVars(combinations, lb = 0.0, vtype=GRB.CONTINUOUS, name="Second_Doses")
    stocks = m.addVars(combinations, lb=0.0, vtype=GRB.CONTINUOUS, name="Stocks")

    m.addConstrs( (first_doses.sum('*',j) + second_doses.sum('*', j) <= Capacity[j] for j in range(0, T)))  

    m.addConstrs( (first_doses[j, i] == second_doses[j,i+Delta] for j, i in combinations if i < T-Delta ))
    m.addConstrs( (first_doses[j, i] == 0 for j, i in combinations if i >= T-Delta ))

    m.addConstrs( second_doses[j,i] == 0 for j, i in combinations if i < Delta)
    
    m.addConstrs ( (first_doses[j,i]  + 0 + stocks[j,i] == B[j,i] + 0 for j, i in combinations if i == 0))
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i] == B[j,i] + stocks[j,i-1] for j, i in combinations if i >= 1 and i < Delta))
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta] + stocks[j,i] == B[j,i] + stocks[j, i-1] for j, i in combinations if i >= Delta))

    m.addConstrs( (stocks[j,i] == 0 for j, i in combinations if i == T-1 ) )

    m.setObjective((second_doses.prod(time_frame)), GRB.MINIMIZE)

    print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        
        print("\n***** Verbose solution printing *****\n")
        
        m.printAttr("X")

        print("\n***** Solution's values list printing *****\n")
        for i in range(0, len(original_B)):

            first_doses_values = resultList[T*i:T*(i+1)]
            print("First_doses values_" + list(original_B)[i] + ": " + str(first_doses_values))

        for i in range(0, len(original_B)):

            second_doses_values = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
            print("Second_doses values_" + list(original_B)[i] + ": " + str(second_doses_values))
        
        for i in range(0, len(original_B)):

            stock_values = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]
            print("Stocks values_" + list(original_B)[i] + ": " + str(stock_values))

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
    
    # Dictionary for the time slot, needed in the objective function

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

def euristic_with_stocks(T, B, Delta):

    print("\n--------------------\n")
    print("***** Test euristic with stocks *****\n")

    second_doses_values = [0.0] * T
    object_function_value = 0

    if (sum(B[:len(B)-Delta]) > int(sum(B)/2)):

        second_doses_values[T-1] = int(sum(B)/2)
        
        for j in range(0, len(second_doses_values)):
            object_function_value += second_doses_values[j] * (j+ Delta + 1)

        print("Second_doses values: " + str(second_doses_values))
        print("Object function value: " + str(object_function_value))
        print("Unused doses: " + str(sum(B) - sum(second_doses_values)*2))

    else:
        print("No solutions available")

'''
def euristic_with_stocks(T, B, Delta):

    print("\n--------------------\n")
    print("***** Test euristic with stocks - Draft *****\n")

    second_doses_values = []
    object_function_value = 0
    stocks = 0

    #       { 0, i < Delta
    # y_i = 
    #       { min (b_{i-Delta}, b_i), i >= Delta

    for i in range(Delta, len(B)):
        
        if i == len(B)-1:
            print(B)
            second_doses_values.append((B[i-Delta]+B[i])/2)
        
        else:

            stocks = abs(B[i-Delta]-B[i])
            min_value = min(B[i-Delta], B[i])

            second_doses_values.append(min_value)
            if B[i-Delta] - stocks >= 0:
                B[i-Delta] -= stocks
                B[i-Delta+1] += B[i-Delta]

            B[i-Delta] -= min_value
            B[i] -= min_value


    for j in range(0, len(second_doses_values)):
        object_function_value += second_doses_values[j] * (j+ Delta + 1)
    
    print("Second_doses values: " + str(second_doses_values))

    print("Object function value: " + str(object_function_value))


    print("Unused doses: " + str(B_value - sum(second_doses_values)*2))
'''

if __name__ == "__main__":
    optimize_test_basic(5, [10,6,8,4,2], 1)
    #optimize_test_capacity(5, [40,30,20,10,20], 3, [100,100,100,100,100])
    #optimize_test_capacity_multiple_vaccines(5, {'Astrazeneca': [40,30,20,10,20], 'Pfizer':[40,30,20,10,20] }, 3, [100,100,100,100,100])
    #euristic(5, [40,30,20,10,20], 3)
    euristic_with_stocks(5, [10,6,8,4,2], 1)