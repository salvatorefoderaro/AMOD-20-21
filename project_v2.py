import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os

PRINT = False
DELTA = {'Pfizer': 21, 'Moderna': 28, 'Astrazeneca': 40}
CSV_INPUT_FOLDER = "csv"
CSV_OUTPUT_FOLDER = "csv_solution_v2"

def optimize_test_capacity_multiple_vaccines(T, B, Delta, Capacity):

    '''
    B = {"Vaccine name": [B_List], "Vaccine name": [B_List] }
    '''

    # print("\n***** Gurobipy solver - Capacity Multiple Vaccines*****\n")

    m = gp.Model("vaccinations")

    if PRINT == False:
        m.Params.LogToConsole = 0

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

    m.addConstrs( (first_doses[j, i] == second_doses[j,i+Delta[j]] for j, i in combinations if i < T-Delta[j] ))
    m.addConstrs( (first_doses[j, i] == 0 for j, i in combinations if i >= T-Delta[j] ))

    m.addConstrs( second_doses[j,i] == 0 for j, i in combinations if i < Delta[j])
    
    m.addConstrs( (first_doses[j,i]  + 0 + stocks[j,i] == B[j,i] + 0 for j, i in combinations if i == 0))
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i] == B[j,i] + stocks[j,i-1] for j, i in combinations if i >= 1 and i < Delta[j]))
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta[j]] + stocks[j,i] == B[j,i] + stocks[j, i-1] for j, i in combinations if i >= Delta[j]))

    m.addConstrs( (stocks[j,i] == 0 for j, i in combinations if i == T-1 ) )

    m.setObjective((second_doses.prod(time_frame)), GRB.MINIMIZE)

    # print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        
        # print("\n***** Verbose solution printing *****\n")
        
        m.printAttr("X")

        # print("\n***** Solution's values list printing *****\n")

        second_doses_dict = {}
        stocks_dict = {}

        for i in range(0, len(original_B)):

            first_doses_values = resultList[T*i:T*(i+1)]
            #print("First_doses values_" + list(original_B)[i] + ": " + str(first_doses_values))

        for i in range(0, len(original_B)):

            second_doses_dict[list(original_B)[i]] = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
            second_doses_values = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
            #print("Second_doses values_" + list(original_B)[i] + ": " + str(second_doses_values))
        
        for i in range(0, len(original_B)):
            stocks_dict[list(original_B)[i]] = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]
            stock_values = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]
            #print("Stocks values_" + list(original_B)[i] + ": " + str(stock_values))

        object_function_value = m.objVal
        # print("Object_function value: " + str(object_function_value))
        return[second_doses_dict, stocks_dict, object_function_value]
    else:
        print("\n***** No solutions found *****")

def get_column_from_df(df, column_name):
    return df[column_name].values.tolist()

def add_column_to_df(df, values_list, new_column_name):
    df[new_column_name] = values_list
    return df

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
    
    optimal_result = [[], []]
    heuristic_result = []
    result_difference = []

    file_list = os.listdir(CSV_INPUT_FOLDER)

    for i in range(0, len(os.listdir(CSV_INPUT_FOLDER)) -2 ):

        data = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i],  index_col=0) 
        data_1 = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i+1],  index_col=0) 
        data_2 = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i+2],  index_col=0) 

        b_list_0 = get_column_from_df(data, "ndosi")
        b_list_1 = get_column_from_df(data_1, "ndosi")
        b_list_2 = get_column_from_df(data_2, "ndosi")

        b_list = {'Pfizer': b_list_0, 'Moderna': b_list_1, 'Astrazeneca': b_list_2}
        total_capacity = sum(b_list_0) + sum(b_list_1) + sum(b_list_2) 
        
        capacity = []
        capacity.append([5 * int(total_capacity/len(b_list))]*len(b_list_0))
        capacity.append([10 * int(total_capacity/len(b_list))]*len(b_list_0))

        index = 0
        for k in capacity:
        # result = optimal_solution(len(b_list), b_list, DELTA)
            result = optimize_test_capacity_multiple_vaccines(len(b_list_0), b_list, DELTA, k)
            if index == 0:
                name = "5c"
            else:
                name = "10c"

            for j in result[0]:
                data = add_column_to_df(data, result[0][j], name + "_second_doses_" + j)
                data = add_column_to_df(data, result[1][j], name + "_stock_values_" + j)

            optimal_result[index].append("{:.2f}".format(result[2]))
            index = index + 1
        #heu_result = heuristic(len(b_list), b_list, DELTA)
        #heuristic_result.append("{:.2f}".format(heu_result))
        #result_difference.append("{:.2f}".format((heu_result-result[2])/(result[2]) *100))
        
        data.to_csv(CSV_OUTPUT_FOLDER + "/solution_" + file_list[i] )
    
    instances = np.arange(1, len(optimal_result[0])+1, 1).tolist()
    df = pd.DataFrame(instances, columns= ['instance'])
    df['5c_optimal_value'] = optimal_result[0]
    df['10c_optimal_value'] = optimal_result[1]
    #df['heuristic_value'] = heuristic_result
    #df['result_difference_%'] = result_difference
    df.to_csv("comparison_result_optimal_heuristic_v2.csv", index=0)