import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os
from statistics import mean
import random

PRINT = False
DELTA = 21
CSV_INPUT_FOLDER = "input_csv"
CSV_OUTPUT_FOLDER = "csv_solution_v1"
STEP = 5

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
                
        return [second_doses_values, stock_values, object_function_value]

    else:
        if PRINT == True:
            print("\n***** No solutions found *****")

def heuristic(T, B, Delta):
    
    x = [0.01]*(T+1)
    y = [0.0]*(T+1)
    s = [0.0] + B[:]
    B = [0.0] + B

    object_function_value = 0

    x[T-Delta] = int(sum(B)/2)
    y[T] = x[T-Delta]
    
    s[T] = 0
    s[T-1] = int(sum(B)/2) - B[T]
    s[T-Delta] = int(sum(B[:T-Delta+1])) - x[T-Delta]

    t = T-Delta

    while(t > 2):
        lista = []
        if x[t] != 0.01:
            lista.append(x[t])
        lista.append(s[t+Delta-1])
        lista.append(B[t-1])

        c = min(lista)
        
        x[t-1] = int(x[t-1]) + c
        y[t+Delta-1] = int(y[t+Delta-1]) + c
        
        x[t] = int(x[t]) - c
        y[t+Delta] = int(y[t+Delta]) - c
        s[t-1] = int(s[t-1]) - c
        s[t+Delta-1] = int(s[t+Delta-1]) -c
       
        t -= 1

    for j in range(0, len(y)):
        object_function_value += y[j] * (j)

    return [object_function_value, y]
    
if __name__ == "__main__":
    
    optimal_result = []
    heuristic_result = []
    result_difference = []
    instances = []
    avg_heuristic = []
    avg_optimal = []

    for i in os.listdir(CSV_INPUT_FOLDER):

        print("Processing instance: " + i)
        
        data = pd.read_csv(CSV_INPUT_FOLDER + "/" + i,  index_col=0) 
        b_list = get_column_from_df(data, "ndosi")
        result = optimal_solution(len(b_list), b_list, DELTA)

        data = add_column_to_df(data, result[0], "second_doses")
        data = add_column_to_df(data, result[1], "stock_values")
        optimal_result.append(round( result[2]/(sum(b_list)/2), 4))
        
        heu_result = heuristic(len(b_list), b_list, DELTA)
        heuristic_result.append(round(heu_result[0]/(sum(b_list)/2), 4))
        result_difference.append(round((heu_result[0]-result[2])/(result[2]) *100, 4))
        
        data.to_csv(CSV_OUTPUT_FOLDER + "/solution_" + i )
    
    for j in range(0, 1000+1-STEP, STEP):
        instances.append(j+STEP)
        avg_heuristic.append(round( mean(heuristic_result[j:j+STEP]), 4))
        avg_optimal.append(round( mean(optimal_result[j:j+STEP]), 4))

    df = pd.DataFrame(instances, columns= ['instance'])
    df['optimal_value'] = avg_optimal
    df['heuristic_value'] = avg_heuristic

    df.to_csv("result_v1.csv", index=0)