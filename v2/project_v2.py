import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os

PRINT = False
DELTA = {'Pfizer': 21, 'Moderna': 28, 'Astrazeneca': 78}
CSV_INPUT_FOLDER = "input_csv"
CSV_OUTPUT_FOLDER = "csv_solution_v2"
T = 180

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

    first_doses = m.addVars(combinations, lb=0.0, vtype=GRB.INTEGER, name="First_Doses")
    second_doses = m.addVars(combinations, lb = 0.0, vtype=GRB.INTEGER, name="Second_Doses")
    stocks = m.addVars(combinations, lb=0.0, vtype=GRB.INTEGER, name="Stocks")

    m.addConstrs( (first_doses.sum('*',j) + second_doses.sum('*', j) <= Capacity for j in range(0, T)))  

    m.addConstrs( (first_doses[j, i] == second_doses[j,i+Delta[j]] for j, i in combinations if i < T-Delta[j] ))
    m.addConstrs( (first_doses[j, i] == 0 for j, i in combinations if i >= T-Delta[j] ))

    m.addConstrs( second_doses[j,i] == 0 for j, i in combinations if i < Delta[j])
    
    m.addConstrs( (first_doses[j,i]  + 0 + stocks[j,i] == B[j,i] + 0 for j, i in combinations if i == 0))
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i] == B[j,i] + stocks[j,i-1] for j, i in combinations if i >= 1 and i < Delta[j]))
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta[j]] + stocks[j,i] == B[j,i] + stocks[j, i-1] for j, i in combinations if i >= Delta[j]))

    # m.addConstrs( (stocks[j,i] == 0 for j, i in combinations if i == T-1 ) )

    m.setObjective((second_doses.prod(time_frame) + 1000 * ( stocks['Pfizer', T-1] + stocks['Moderna', T-1] + stocks['Astrazeneca', T-1])), GRB.MINIMIZE)

    # print ("\n\n***** Optimize log *****\n\n")

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
        
        # print("\n***** Verbose solution printing *****\n")
        
        m.printAttr("X")

        # print("\n***** Solution's values list printing *****\n")

        first_doses_dict = {}
        second_doses_dict = {}
        stocks_dict = {}

        for i in range(0, len(original_B)):
            
            first_doses_dict[list(original_B)[i]] = resultList[T*i:T*(i+1)]
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
        return[first_doses_dict, second_doses_dict, stocks_dict, object_function_value]
    else:
        print("\n***** No solutions found *****")

def get_column_from_df(df, column_name):
    return df[column_name].values.tolist()

def add_column_to_df(df, values_list, new_column_name):
    df[new_column_name] = values_list
    return df

if __name__ == "__main__":
    
    optimal_result = {"c": [], "2c": [], "3c": [],"4c": [], "5c": [], "6c": [], "7c": [], "8c": [], "9c": [], "10c": [], "11c": []}
    heuristic_result = []
    result_difference = []

    file_list = os.listdir(CSV_INPUT_FOLDER)

    instances = np.arange(1, T + 1, 1).tolist()

    for i in range(0, len(os.listdir(CSV_INPUT_FOLDER)) -2-900 ):

        print("Processing instance: " + str(i))

        df = pd.DataFrame(instances, columns= ['instance'])

        # Read the csv files
        data = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i],  index_col=0) 
        data_1 = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i+1],  index_col=0) 
        data_2 = pd.read_csv(CSV_INPUT_FOLDER + "/" + file_list[i+2],  index_col=0) 

        # Get column from the csv files
        b_list_0 = get_column_from_df(data, "ndosi")
        b_list_1 = get_column_from_df(data_1, "ndosi")
        b_list_2 = get_column_from_df(data_2, "ndosi")

        df["arrival_pfizer"] = b_list_0
        df["arrival_moderna"] = b_list_1
        df["arrival_astrazeneca"] = b_list_2

        # Dictionary for multiple vaccines
        b_list = {'Pfizer': b_list_0, 'Moderna': b_list_1, 'Astrazeneca': b_list_2}
        
        # List of different capacity
        total_capacity = sum(b_list_0) + sum(b_list_1) + sum(b_list_2)

        capacity = {}
        capacity["c"] =   1 * int(total_capacity/180)
        capacity["2c"] =  1.1 * int(total_capacity/180) 
        capacity["3c"] =  1.2 * int(total_capacity/180) 
        capacity["4c"] =  1.3 * int(total_capacity/180) 
        capacity["5c"] =  1.4 * int(total_capacity/180) 
        capacity["6c"] =  1.5 * int(total_capacity/180)
        capacity["7c"] =  1.6 * int(total_capacity/180) 
        capacity["8c"] =  1.7 * int(total_capacity/180) 
        capacity["9c"] =  1.8 * int(total_capacity/180) 
        capacity["10c"] = 1.9 * int(total_capacity/180) 
        capacity["11c"] = 2.0 * int(total_capacity/180) 

        # Calculate the optimal value

        for u in capacity:

            result = optimize_test_capacity_multiple_vaccines(len(b_list_0), b_list, DELTA, capacity[u])

            for j in result[0]:
                df["first_doses_" + j] = result[0][j]
                df["second_doses_" + j] = result[1][j]
                df["stock_values_" + j] = result[2][j]

            optimal_result[u].append(result[3]/180)
            df["capacity"] = [capacity[u]] * 180
        
            # Write the solution to the single file
            # df.to_csv(CSV_OUTPUT_FOLDER + "/capacity_" + u +"/solution_" + file_list[i] )

    instances = np.arange(1, len(optimal_result["c"]) + 1-900, 1).tolist()

    # Write the solution to the summary file
    df = pd.DataFrame(instances, columns= ['instance'])
    df['c_optimal_value'] = optimal_result["c"]
    df['1.1c_optimal_value'] = optimal_result["2c"]
    df['1.2c_optimal_value'] = optimal_result["3c"]
    df['1.3c_optimal_value'] = optimal_result["4c"]
    df['1.4c_optimal_value'] = optimal_result["5c"]
    df['1.5c_optimal_value'] = optimal_result["6c"]
    df['1.6c_optimal_value'] = optimal_result["7c"]
    df['1.7c_optimal_value'] = optimal_result["8c"]   
    df['1.8c_optimal_value'] = optimal_result["9c"]
    df['1.9c_optimal_value'] = optimal_result["10c"] 
    df['2c_optimal_value'] = optimal_result["11c"] 

    df.to_csv("result_v2.csv", index=0)