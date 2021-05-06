import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os

PRINT = False
DELTA = 21
CSV_INPUT_FOLDER = "csv"
CSV_OUTPUT_FOLDER = "csv_solution"

def get_column_from_df(df, column_name):
    return df[column_name].values.tolist()

def add_column_to_df(df, values_list, new_column_name):
    df[new_column_name] = values_list
    return df

def optimal_solution(T, B, Delta):

    if PRINT == True:
        print("\n***** Gurobipy solver - Basic*****\n")

    # Declare the model
    
    m = gp.Model("vaccinations")

    if PRINT == False:
        m.Params.LogToConsole = 0

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

    # print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        first_doses_values = resultList[:T]
        second_doses_values = resultList[T:2*T]
        stock_values = resultList[2*T:]
        object_function_value = m.objVal
        
        if PRINT == True:
            m.printAttr("X")
            print("\n***** Verbose solution printing *****\n")
            print("\n***** Solution's values list printing *****\n")
            print("First_doses values: " + str(first_doses_values))
            print("Second_doses values: " + str(second_doses_values))
            print("Stocks values: " + str(stock_values))
            print("Object_function value: " + str(object_function_value))
    
        return [second_doses_values, stock_values, object_function_value/T]

    else:
        if PRINT == True:
            print("\n***** No solutions found *****")

def heuristic(T, B, Delta):

    x = [0.0]*T
    y = [0.0]*T
    s = B[:]

    object_function_value = 0

    x[T-Delta-1] = int(sum(B)/2)
    y[T-1] = x[T-Delta-1]
    
    s[T-1] = 0
    s[T-Delta-1] = sum(B[:T-Delta]) - x[T-Delta-1]

    t = T-Delta-1

    while(t > 0):

        c = min(s[t+Delta-1], x[t], B[t-1])
        
        x[t-1] += c
        y[t+Delta-1] += c
        
        x[t] -= c
        y[t+Delta] -= c
        s[t-1] -= c
        s[t+Delta-1] -= c
       
        t -= 1

    for j in range(0, len(y)):
        object_function_value += y[j] * (j + 1)
    
    return object_function_value/T
    
if __name__ == "__main__":
    
    optimal_result = []
    heuristic_result = []
    result_difference = []

    for i in os.listdir(CSV_INPUT_FOLDER):
        
        data = pd.read_csv(CSV_INPUT_FOLDER + "/" + i,  index_col=0) 
        b_list = get_column_from_df(data, "ndosi")
        result = optimal_solution(len(b_list), b_list, DELTA)

        data = add_column_to_df(data, result[0], "second_doses")
        data = add_column_to_df(data, result[1], "stock_values")
        optimal_result.append("{:.2f}".format(result[2]))
        
        heu_result = heuristic(len(b_list), b_list, DELTA)
        heuristic_result.append("{:.2f}".format(heu_result))
        result_difference.append("{:.2f}".format((heu_result-result[2])/(result[2]) *100))
        
        data.to_csv(CSV_OUTPUT_FOLDER + "/solution_" + i )
    
    instances = np.arange(1, len(optimal_result)+1, 1).tolist()
    df = pd.DataFrame(instances, columns= ['instance'])
    df['optimal_value'] = optimal_result
    df['heuristic_value'] = heuristic_result
    df['result_difference_%'] = result_difference
    df.to_csv("comparison_result_optimal_heuristic.csv", index=0)