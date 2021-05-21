import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os
from pathlib import Path
from statistics import mean
import sys
import shutil

PRINT = False
DELTA = {'Pfizer': 21, 'Moderna': 28, 'Astrazeneca': 78}
CSV_INPUT_FOLDER = "input_csv"
CSV_OUTPUT_FOLDER = "csv_solution_v2"
T = 180
NUMBER_OF_ELEMENT = 30
LIMIT_CAPACITY = 8
INCREMENT = 0.5
LAST_DAY_VACCINES = True

def optimize_test_capacity_multiple_vaccines_robust(T, Delta, Capacity, Instance):

    m = gp.Model("vaccinations")

    Capacity = 2 * Capacity

    m.Params.LogToConsole = 0
    vaccine = {"Pfizer", "Moderna", "Astrazeneca"}
    dicti = {}
    dict_B = {}
    dict_stocks = {}

    for i in vaccine:
        for j in range(0, T):
            dicti[(i, j)] = [j+1]
            for k in Instance:
                dict_stocks[(i,j,k)] = [j+1]

    combinations, time_frame = gp.multidict(dicti)
    combinations_stock, time_frame = gp.multidict(dict_stocks)

    first_doses = m.addVars(combinations, lb=0.0, vtype=GRB.INTEGER, name="First_Doses")
    second_doses = m.addVars(combinations, lb = 0.0, vtype=GRB.INTEGER, name="Second_Doses")
    stocks = m.addVars(combinations_stock, lb=0.0, vtype=GRB.INTEGER, name="Stocks")
    z = m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name="Z")

    m.addConstrs( (first_doses.sum('*',j) + second_doses.sum('*', j) <= Capacity for j in range(0, T)))  

    m.addConstrs( (first_doses[j, i] == second_doses[j,i+Delta[j]] for j, i in combinations if i < T-Delta[j] ))
    m.addConstrs( (first_doses[j, i] == 0 for j, i in combinations if i >= T-Delta[j] ))

    m.addConstrs( second_doses[j,i] == 0 for j, i in combinations if i < Delta[j])
    #m.addConstrs( stocks[j, i, k] <= 10000  for j, i, k in combinations_stock if i == T-1)
    
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i,k] == Instance[k][j][i] + 0 for j, i, k in combinations_stock if i == 0))
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i,k] == Instance[k][j][i] + stocks[j,i-1,k] for j, i, k in combinations_stock if i >= 1 and i < Delta[j]))
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta[j]] + stocks[j,i,k] == Instance[k][j][i] + stocks[j,i-1,k] for j, i, k in combinations_stock if i >= 1 and i >= Delta[j]))

    for k in Instance:
        m.addConstr( z >= (( 2 / (Instance[k]["TotalB"]) * (gp.quicksum(second_doses[j,i] * i for j,i in combinations) + ( Instance[k]["stock_values_Pfizer"]*(180 + 21) + (Instance[k]["stock_values_Moderna"])*(180 + 28) + (Instance[k]["stock_values_Astrazeneca"])*(180 + 78)))) - Instance[k]["Result"] ) )

    m.setObjective(  z , GRB.MINIMIZE)
    # m.addConstrs( stocks[j, i, k] <= 15000 for j, i, k in stocks_combinations if i == T-1)

    #m.addConstrs ( z >=  ( ( (gp.quicksum(second_doses[j,i] * i   for j,i in combinations))) )   for k in Instance )

    #m.setObjective(  gp.quicksum(second_doses[j,i] * i  + Instance[0]['stock_values_Astrazeneca'] + Instance[0]['stock_values_Moderna'] + Instance[0]['stock_values_Pfizer'] for j,i in combinations ), GRB.MINIMIZE)

    m.optimize()

    if (m.solCount > 0):

        print(m.objVal)

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
            
        first_doses_dict = {}
        second_doses_dict = {}
        stocks_dict = {}

        '''for i in range(0, len(original_B)):
            
            first_doses_dict[list(original_B)[i]] = resultList[T*i:T*(i+1)]
            first_doses_values = resultList[T*i:T*(i+1)]

        for i in range(0, len(original_B)):

            second_doses_dict[list(original_B)[i]] = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
            second_doses_values = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
        
        for i in range(0, len(original_B)):
            stocks_dict[list(original_B)[i]] = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]
            stock_values = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]'''

        object_function_value = m.objVal
        return[first_doses_dict, second_doses_dict, stocks_dict, object_function_value]
    else:
        print("\n***** No solutions found *****")

def optimize_test_capacity_multiple_vaccines(T, B, Delta, Capacity):

    m = gp.Model("vaccinations")

    #if PRINT == False:

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
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta[j]] + stocks[j,i] == B[j,i] + stocks[j, i-1] for j, i in combinations if i >= 1 and i >= Delta[j]))

    m.setObjective(  (gp.quicksum(second_doses[j,i] * i for j,i in combinations) + 1000 * ( stocks['Pfizer', T-1] + stocks['Moderna', T-1] + stocks['Astrazeneca', T-1])), GRB.MINIMIZE)

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
            
        first_doses_dict = {}
        second_doses_dict = {}
        stocks_dict = {}

        for i in range(0, len(original_B)):
            
            first_doses_dict[list(original_B)[i]] = resultList[T*i:T*(i+1)]
            first_doses_values = resultList[T*i:T*(i+1)]

        for i in range(0, len(original_B)):

            second_doses_dict[list(original_B)[i]] = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
            second_doses_values = resultList[T*(i+len(original_B)):T*(i+1+len(original_B))]
        
        for i in range(0, len(original_B)):
            stocks_dict[list(original_B)[i]] = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]
            stock_values = resultList[T*(i+2*len(original_B)):T*(i+1+2*len(original_B))]

        object_function_value = m.objVal
        return[first_doses_dict, second_doses_dict, stocks_dict, object_function_value]
    else:
        print("\n***** No solutions found *****")

def get_column_from_df(df, column_name):
    return df[column_name].values.tolist()

def add_column_to_df(df, values_list, new_column_name):
    df[new_column_name] = values_list
    return df

def heuristic_v2_sum(first_doses_list, delta, t, q):
    value_sum = 0

    for i in range(1, q+1):
        if t + i - delta >= 0 and t + i - delta <= 180:
            value_sum += first_doses_list[t + i - delta]
    return value_sum

def heuristic_v2_risk(b_list, capacity, q):

    vaccines = { "Pfizer": [21, [0]*180, [0]*180], "Moderna": [28, [0]*180, [0]*180], "Astrazeneca": [78, [0]*180, [0]*180] }
    total_arrival = {"Pfizer": sum(b_list["Pfizer"]), "Moderna": sum(b_list["Moderna"]), "Astrazeneca": sum(b_list["Astrazeneca"])   }
    doses_somministrated = {"Pfizer": 0, "Moderna": 0, "Astrazeneca": 0}
    arrival_sum = total_arrival["Pfizer"] + total_arrival["Moderna"] + total_arrival["Astrazeneca"]

    sum_total_arrival = total_arrival["Pfizer"] + total_arrival["Moderna"] + total_arrival["Astrazeneca"]
    object_function = 0 
    total_second_doses = 0
    index = 0
    sum_penality = 0
    count_negative = 0

    negative_stocks = 0
    total_stocks = 0
    utilization = 0
    min_negative_stocks = 0

    for t in range(0, 180):
        index = 0
        remain_capacity = capacity
        for u in vaccines:

            first_doses_somministrated = 0
            second_doses_somministrated = 0
            delta = vaccines[u][0]
            first_doses = vaccines[u][1]
            stocks = vaccines[u][2]
            arrival = b_list[u]

            if t + delta < 180:

                if t == 0:
                    first_doses_somministrated = arrival[t]
                elif t-delta < 0:
                    first_doses_somministrated = stocks[t-1] + arrival[t]
                else:
                    first_doses_somministrated = stocks[t-1] - first_doses[t-delta] - heuristic_v2_sum(first_doses, delta, t, q) + arrival[t]

                if first_doses_somministrated < 0:
                    first_doses_somministrated = 0

                if t-delta >= 0:
                    second_doses_somministrated = first_doses[t-delta]
                else:
                    second_doses_somministrated = 0

                first_doses_somministrated = min(first_doses_somministrated , (capacity * 0.33) - second_doses_somministrated)

                object_function += (t+1+delta)*first_doses_somministrated
                utilization += first_doses_somministrated
                doses_somministrated[u] += first_doses_somministrated
                vaccines[u][1][t] = first_doses_somministrated

                if t == 0:
                    vaccines[u][2][t] = - first_doses_somministrated  + arrival[t]
                elif t - delta < 0:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated + arrival[t]
                else:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated - first_doses[t-delta] + arrival[t]
            else:
                vaccines[u][2][t] = stocks[t-1]  - first_doses[t-delta]  + arrival[t]

            if vaccines[u][2][t] < 0:
                negative_stocks += vaccines[u][2][t]
                count_negative += 1
                min_negative_stocks = min(min_negative_stocks, vaccines[u][2][t])
            '''else:
                utilization += first_doses[t-delta]
                object_function += (t+1) * first_doses[t-delta]'''
            total_stocks = vaccines[u][2][t]

    total_dose_somministrated =  0
    total_dose_somministrated += doses_somministrated["Pfizer"] * 2
    total_dose_somministrated += doses_somministrated["Astrazeneca"] * 2
    total_dose_somministrated += doses_somministrated["Moderna"] * 2

    stocks_at_end = 0
    stocks_at_end += total_arrival["Pfizer"] - doses_somministrated["Pfizer"] * 2
    stocks_at_end += total_arrival["Moderna"] - doses_somministrated["Moderna"] * 2
    stocks_at_end += total_arrival["Astrazeneca"] - doses_somministrated["Astrazeneca"] * 2

    if LAST_DAY_VACCINES:

        object_function += (total_arrival["Pfizer"] - doses_somministrated["Pfizer"] * 2) * (180+21)
        object_function += (total_arrival["Moderna"] - doses_somministrated["Moderna"] * 2) * (180+28)
        object_function += (total_arrival["Astrazeneca"] - doses_somministrated["Astrazeneca"] * 2) * (180+78)
        return [ (object_function ) / ((total_dose_somministrated + stocks_at_end)/2) , 0, abs(negative_stocks)/arrival_sum, count_negative/ (3*180), (utilization)/sum_total_arrival, min_negative_stocks]

    else:
        return [ (object_function ) / (total_dose_somministrated) , stocks_at_end, abs(negative_stocks)/arrival_sum, (count_negative)/ (3*180), (utilization)/sum_total_arrival, min_negative_stocks]


def heuristic_v2(b_list, capacity):

    vaccines = { "Pfizer": [21, [0]*180, [0]*180], "Moderna": [28, [0]*180, [0]*180], "Astrazeneca": [78, [0]*180, [0]*180] }
    stocks = []
    total_arrival = {"Pfizer": sum(b_list["Pfizer"]), "Moderna": sum(b_list["Moderna"]), "Astrazeneca": sum(b_list["Astrazeneca"])   }
    doses_somministrated = {"Pfizer": 0, "Moderna": 0, "Astrazeneca": 0}
    sum_total_arrival = total_arrival["Pfizer"] + total_arrival["Moderna"] + total_arrival["Astrazeneca"]

    object_function = 0 
    total_second_doses = 0
    index = 0
    sum_penality = 0

    for t in range(0, 180):
        
        index = 0
        remain_capacity = capacity

        for u in vaccines:
            delta = vaccines[u][0]
            first_doses = vaccines[u][1]
            stocks = vaccines[u][2]
            arrival = b_list[u]

            if t + delta < 180:

                if t == 0:
                    first_doses_somministrated = arrival[t]
                elif t-delta < 0:
                    first_doses_somministrated = stocks[t-1] + arrival[t]
                else:
                    first_doses_somministrated = stocks[t-1] - first_doses[t-delta] + arrival[t]

                first_doses_somministrated = int(first_doses_somministrated / 2)

                if t-delta >= 0:
                    second_doses_somministrated = first_doses[t-delta]
                else:
                    second_doses_somministrated = 0
    
                first_doses_somministrated = min(first_doses_somministrated, (capacity * 0.33) - second_doses_somministrated)

                total_second_doses += first_doses_somministrated
                object_function += (t+1+delta)*first_doses_somministrated
                vaccines[u][1][t] = first_doses_somministrated
                doses_somministrated[u] += first_doses_somministrated
            
                if t == 0:
                    vaccines[u][2][t] = - first_doses_somministrated  + arrival[t]
                elif t - delta < 0:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated + arrival[t]
                else:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated - first_doses[t-delta] + arrival[t]
            else:
                vaccines[u][2][t] = stocks[t-1]  - first_doses[t-delta]  + arrival[t]

    total_dose_somministrated =  0
    total_dose_somministrated += doses_somministrated["Pfizer"] * 2
    total_dose_somministrated += doses_somministrated["Astrazeneca"] * 2
    total_dose_somministrated += doses_somministrated["Moderna"] * 2

    stocks_at_end = 0
    stocks_at_end += total_arrival["Pfizer"] - doses_somministrated["Pfizer"] * 2
    stocks_at_end += total_arrival["Moderna"] - doses_somministrated["Moderna"] * 2
    stocks_at_end += total_arrival["Astrazeneca"] - doses_somministrated["Astrazeneca"] * 2

    if stocks_at_end - (vaccines["Astrazeneca"][2][179] +  vaccines["Pfizer"][2][179] +  vaccines["Moderna"][2][179]) > 100:
        print("Errore")
        input()

    if LAST_DAY_VACCINES:

        object_function += (total_arrival["Pfizer"] - doses_somministrated["Pfizer"] * 2) * (180+21)
        object_function += (total_arrival["Moderna"] - doses_somministrated["Moderna"] * 2) * (180+28)
        object_function += (total_arrival["Astrazeneca"] - doses_somministrated["Astrazeneca"] * 2) * (180+78)

        return [ (object_function ) / ((total_dose_somministrated + stocks_at_end)/2) , 0, (total_dose_somministrated/2)/sum_total_arrival]
    
    else:
        return [ (object_function ) / (total_dose_somministrated) , stocks_at_end, (total_dose_somministrated/2)/sum_total_arrival]

if __name__ == "__main__":

    shutil.rmtree("csv_solution_v2")
    os.mkdir("csv_solution_v2")

    if len(sys.argv) > 1:
        if sys.argv[1] == "True":
            LAST_DAY_VACCINES = True
            INCREMENT = float(sys.argv[2])
            NUMBER_OF_ELEMENT = int(sys.argv[3])
        else:
            LAST_DAY_VACCINES = False
            INCREMENT = float(sys.argv[2])
            NUMBER_OF_ELEMENT = int(sys.argv[3])
    
    robust_optimization = {}
    penality_optimal_result = {}
    optimal_result = {}

    min_negative_stocks_risk = {}
    min_negative_stocks_risk_7 = {}
    min_negative_stocks_risk_14 = {}
    min_negative_stocks_risk_21 = {}

    feasible_s_pfizer = 0
    feasible_s_moderna = 0
    feasible_s_astrazeneca = 0
    
    heuristic_result = {}
    heuristic_risk_result = {}
    heuristic_risk_result_7 = {}
    heuristic_risk_result_14 = {}
    heuristic_risk_result_21 = {}

    utilization_T = {}
    utilization_T_heuristic = {}
    utilization_T_heuristic_risk = {}
    utilization_T_heuristic_risk_7 = {}
    utilization_T_heuristic_risk_14 = {}
    utilization_T_heuristic_risk_21 = {}

    penality_heuristic = {}
    penality_heuristic_risk = {}
    penality_heuristic_risk_7 = {}
    penality_heuristic_risk_14 = {}
    penality_heuristic_risk_21 = {}

    heuristic_risk_count_negative = {}
    heuristic_risk_count_negative_7 = {}
    heuristic_risk_count_negative_14 = {}
    heuristic_risk_count_negative_21 = {}

    heuristic_risk_negative_arrival = {}
    heuristic_risk_negative_arrival_7 = {}
    heuristic_risk_negative_arrival_14 = {}
    heuristic_risk_negative_arrival_21 = {}

    num = 1
    while(num < (LIMIT_CAPACITY)):
        penality_optimal_result[str(num) + " c"] = []
        optimal_result[str(num) + " c"] = []

        min_negative_stocks_risk[str(num) + " c"] = 0
        min_negative_stocks_risk_7[str(num) + " c"] = 0
        min_negative_stocks_risk_14[str(num) + " c"] = 0
        min_negative_stocks_risk_21[str(num) + " c"] = 0
        
        heuristic_result[str(num) + " c"] = []
        heuristic_risk_result[str(num) + " c"] = []
        heuristic_risk_result_7[str(num) + " c"] = []
        heuristic_risk_result_14[str(num) + " c"] = []
        heuristic_risk_result_21[str(num) + " c"] = []

        penality_heuristic[str(num) + " c"] = []
        penality_heuristic_risk[str(num) + " c"] = []
        penality_heuristic_risk_7[str(num) + " c"] = []
        penality_heuristic_risk_14[str(num) + " c"] = []
        penality_heuristic_risk_21[str(num) + " c"] = []

        utilization_T[str(num) + " c"] = []
        utilization_T_heuristic[str(num) + " c"] = []
        utilization_T_heuristic_risk[str(num) + " c"] = []
        utilization_T_heuristic_risk_7[str(num) + " c"] = []
        utilization_T_heuristic_risk_14[str(num) + " c"] = []
        utilization_T_heuristic_risk_21[str(num) + " c"] = []

        heuristic_risk_count_negative[str(num) + " c"] = []
        heuristic_risk_count_negative_7[str(num) + " c"] = []
        heuristic_risk_count_negative_14[str(num) + " c"] = []
        heuristic_risk_count_negative_21[str(num) + " c"] = []

        heuristic_risk_negative_arrival[str(num) + " c"] = []
        heuristic_risk_negative_arrival_7[str(num) + " c"] = []
        heuristic_risk_negative_arrival_14[str(num) + " c"] = []
        heuristic_risk_negative_arrival_21[str(num) + " c"] = []

        num += INCREMENT
        num = round(num, 2)

    file_list = os.listdir(CSV_INPUT_FOLDER)
    instances = np.arange(1, T + 1, 1).tolist()

    for i in range(0, len(os.listdir(CSV_INPUT_FOLDER)) -2  - (1000 - NUMBER_OF_ELEMENT) ):

        if i == 45:
            continue

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

        b_list = {'Pfizer': b_list_0, 'Moderna': b_list_1, 'Astrazeneca': b_list_2}
        
        # Build the dict with different capacity
        total_capacity = sum(b_list_0) + sum(b_list_1) + sum(b_list_2)
        capacity = {}

        num = 1
        while(num < LIMIT_CAPACITY):
            capacity[str(num) + " c"] =  int(num * int(total_capacity/180))
            num += INCREMENT
            num = round(num, 2)

        feasible_s_pfizer = max(feasible_s_pfizer, min(0, sum(b_list_0) - 2*sum(b_list_0[:180-21])))
        feasible_s_moderna = max(feasible_s_moderna, min(0, sum(b_list_1) - 2*sum(b_list_1[:180-28])))
        feasible_s_astrazeneca = max(feasible_s_astrazeneca, min(0, sum(b_list_2) - 2*sum(b_list_2[:180-78])))

        # For each capacity...
        for u in capacity:

            result = optimize_test_capacity_multiple_vaccines(len(b_list_0), b_list, DELTA, capacity[u])
            heu_result = heuristic_v2(b_list, capacity[u])
            heu_result_risk = heuristic_v2_risk(b_list, capacity[u], 1)
            heu_result_risk_7 = heuristic_v2_risk(b_list, capacity[u], 7)
            heu_result_risk_14 = heuristic_v2_risk(b_list, capacity[u], 14)
            heu_result_risk_21 = heuristic_v2_risk(b_list, capacity[u], 21)

            # Variable to get the penality for the optimal solution
            second_doses_sum = 0
            penality_sum = 0
            new_doses = 0

            robust_optimization[i] = {}

            for j in result[0]:
                df["first_doses_" + j] = result[0][j]
                df["second_doses_" + j] = result[1][j]
                df["stock_values_" + j] = result[2][j]
                robust_optimization[i]["stock_values_" + j] = result[2][j][180-1]
                second_doses_sum += sum(result[1][j])
                penality_sum += result[2][j][180-1]
                if LAST_DAY_VACCINES:
                    new_doses += (result[2][j][180-1]/2) * (180 + DELTA[j]) 

            # Build the capacity list to append at the end of the CSV file
            df["capacity"] = [capacity[u]] * 180
            
            # Calculate the optimal result and append to the CSV file
            solution_penality = 1000 * penality_sum
            
            if LAST_DAY_VACCINES:
                optimal_result_without_penality = (result[3] - solution_penality + new_doses) / (second_doses_sum + penality_sum/2)
            else:
                optimal_result_without_penality = (result[3] - solution_penality) / (second_doses_sum)

            optimal_result[u].append( optimal_result_without_penality )

            robust_optimization[i].update({
                "Result" : optimal_result_without_penality,
                "Pfizer":  b_list_0,
                "Moderna": b_list_1,
                "Astrazeneca": b_list_2,
                "TotalB": sum(b_list_0) + sum(b_list_1) + sum(b_list_2)
                })

            if LAST_DAY_VACCINES:
                penality_optimal_result[u].append(0)
            else:
                penality_optimal_result[u].append( penality_sum )

            utilization_T[u].append( second_doses_sum / total_capacity )

            min_negative_stocks_risk[u] = min(min_negative_stocks_risk[u], heu_result_risk[5])
            min_negative_stocks_risk_7[u] = min(min_negative_stocks_risk_7[u], heu_result_risk_7[5])
            min_negative_stocks_risk_14[u] = min(min_negative_stocks_risk_14[u], heu_result_risk_14[5])
            min_negative_stocks_risk_21[u] = min(min_negative_stocks_risk_21[u], heu_result_risk_21[5])

            # Calculate the heuristic result
            heuristic_result[u].append( heu_result[0] )
            penality_heuristic[u].append( heu_result[1] )
            utilization_T_heuristic[u].append( heu_result[2] )

            heuristic_risk_result[u].append( heu_result_risk[0] )
            penality_heuristic_risk[u].append( heu_result_risk[1] )
            heuristic_risk_negative_arrival[u].append( heu_result_risk[2] )
            heuristic_risk_count_negative[u].append( heu_result_risk[3] )
            utilization_T_heuristic_risk[u].append( heu_result_risk[4] )

            heuristic_risk_result_7[u].append( heu_result_risk_7[0] )
            penality_heuristic_risk_7[u].append( heu_result_risk_7[1] )
            heuristic_risk_negative_arrival_7[u].append( heu_result_risk_7[2] )
            heuristic_risk_count_negative_7[u].append( heu_result_risk_7[3] )
            utilization_T_heuristic_risk_7[u].append( heu_result_risk_7[4] )

            heuristic_risk_result_14[u].append( heu_result_risk_14[0] )
            penality_heuristic_risk_14[u].append( heu_result_risk_14[1] )
            heuristic_risk_negative_arrival_14[u].append( heu_result_risk_14[2] )
            heuristic_risk_count_negative_14[u].append( heu_result_risk_14[3] )
            utilization_T_heuristic_risk_14[u].append( heu_result_risk_14[4] )

            heuristic_risk_result_21[u].append( heu_result_risk_21[0] )
            penality_heuristic_risk_21[u].append( heu_result_risk_21[1] )
            heuristic_risk_negative_arrival_21[u].append( heu_result_risk_21[2] )
            heuristic_risk_count_negative_21[u].append( heu_result_risk_21[3] )
            utilization_T_heuristic_risk_21[u].append( heu_result_risk_21[4] )

            # Write the solution to the single file (create direcitory if not exists)
            output_dir = Path("csv_solution_v2/capacity_" + u)
            output_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(CSV_OUTPUT_FOLDER + "/capacity_" + u + "/solution_" + file_list[i] )

    instances = np.arange(1, len(optimal_result["1 c"]) + 1 - (1000 - NUMBER_OF_ELEMENT), 1).tolist()

    print(min_negative_stocks_risk)
    print(min_negative_stocks_risk_7)
    print(min_negative_stocks_risk_14)
    print(min_negative_stocks_risk_21)


    for i in robust_optimization:
        robust_optimization[i]["sPfizer"] = feasible_s_pfizer
        robust_optimization[i]["sModerna"] = feasible_s_moderna
        robust_optimization[i]["sAstrazeneca"] = feasible_s_astrazeneca

    # optimize_test_capacity_multiple_vaccines_robust(180, DELTA, capacity["1 c"],robust_optimization)

    # Write the solution to the summary file
    df = pd.DataFrame(instances, columns= ['instance'])
    df_summary = pd.DataFrame(instances, columns= ['instance'])
    capacity_list = []

    avg_optimal_value = []
    avg_heuristic_value = []
    avg_heuristic_risk_value = []
    avg_heuristic_risk_value_7 = []
    avg_heuristic_risk_value_14 = []
    avg_heuristic_risk_value_21 = []

    avg_result_difference = []
    avg_result_difference_risk = []
    avg_result_difference_risk_7 = []
    avg_result_difference_risk_14 = []
    avg_result_difference_risk_21 = []

    avg_stocks_optimal = []
    avg_stocks_heuristic = []
    avg_stocks_risk = []
    avg_stocks_risk_7 = []
    avg_stocks_risk_14 = []
    avg_stocks_risk_21 = []

    avg_stocks_difference = []
    avg_stocks_difference_risk = []
    avg_stocks_difference_risk_7 = []
    avg_stocks_difference_risk_14 = []
    avg_stocks_difference_risk_21 = []

    avg_utilization_T = []
    avg_utilization_T_heuristic = []
    avg_utilization_T_heuristic_risk = []
    avg_utilization_T_heuristic_risk_7 = []
    avg_utilization_T_heuristic_risk_14 = []
    avg_utilization_T_heuristic_risk_21 = []

    avg_heuristic_risk_negative_arrival = []
    avg_heuristic_risk_negative_arrival_7 = []
    avg_heuristic_risk_negative_arrival_14 = []
    avg_heuristic_risk_negative_arrival_21 =[]

    avg_heuristic_risk_count_negative = []
    avg_heuristic_risk_count_negative_7 = []
    avg_heuristic_risk_count_negative_14 = []
    avg_heuristic_risk_count_negative_21 = []

    total_min_negative_stocks_risk = []
    total_min_negative_stocks_risk_7 = []
    total_min_negative_stocks_risk_14 = []
    total_min_negative_stocks_risk_21 = []


    for k in optimal_result:
        df[k + " - optimal solution"] = optimal_result[k]
        df[k + " - optimal solution remain stocks"] = penality_optimal_result[k]

        df[k + " - conservative heuristic solution"] = heuristic_result[k]
        df[k + " - conservative heuristic remain stocks"] = penality_heuristic[k]

        df[k + " - heuristic q-1 risk solution"] = heuristic_risk_result[k]
        df[k + " - heuristic q-1 risk remain stocks"] = penality_heuristic_risk[k]

        df[k + " - heuristic q-7 risk solution"] = heuristic_risk_result_7[k]
        df[k + " - heuristic q-7risk remain stocks"] = penality_heuristic_risk_7[k]

        df[k + " - heuristic q14 risk solution"] = heuristic_risk_result_14[k]
        df[k + " - heuristic q-14 risk remain stocks"] = penality_heuristic_risk_14[k]

        df[k + " - heuristic q-21 risk solution"] = heuristic_risk_result_21[k]
        df[k + " - heuristic q-21 risk remain stocks"] = penality_heuristic_risk_21[k]
        
        capacity_list.append(k)

        avg_optimal_value.append(round( mean(optimal_result[k] ), 2))
        avg_heuristic_value.append(round( mean(heuristic_result[k] ), 2))
        avg_heuristic_risk_value.append(round( mean(heuristic_risk_result[k] ), 2))
        avg_heuristic_risk_value_7.append(round( mean(heuristic_risk_result_7[k] ), 2))
        avg_heuristic_risk_value_14.append(round( mean(heuristic_risk_result_14[k] ), 2))
        avg_heuristic_risk_value_21.append(round( mean(heuristic_risk_result_21[k] ), 2))

        avg_stocks_optimal.append(round( mean(penality_optimal_result[k] ), 2))
        avg_stocks_heuristic.append(round( mean(penality_heuristic[k] ), 2))
        avg_stocks_risk.append(round( mean(penality_heuristic_risk[k] ), 2))
        avg_stocks_risk_7.append(round( mean(penality_heuristic_risk_7[k] ), 2))
        avg_stocks_risk_14.append(round( mean(penality_heuristic_risk_14[k] ), 2))
        avg_stocks_risk_21.append(round( mean(penality_heuristic_risk_21[k] ), 2))

        avg_utilization_T.append( round( mean(utilization_T[k] ), 5) )
        avg_utilization_T_heuristic.append( round( mean(utilization_T_heuristic[k] ), 5) )
        avg_utilization_T_heuristic_risk.append( round( mean(utilization_T_heuristic_risk[k] ), 5) )
        avg_utilization_T_heuristic_risk_7.append( round( mean(utilization_T_heuristic_risk_7[k] ), 5) )
        avg_utilization_T_heuristic_risk_14.append( round( mean(utilization_T_heuristic_risk_14[k] ), 5) )
        avg_utilization_T_heuristic_risk_21.append( round( mean(utilization_T_heuristic_risk_21[k] ), 5) )

        avg_heuristic_risk_count_negative.append(round (mean(heuristic_risk_count_negative[k]), 4))
        avg_heuristic_risk_count_negative_7.append(round (mean(heuristic_risk_count_negative_7[k]), 4))
        avg_heuristic_risk_count_negative_14.append(round (mean(heuristic_risk_count_negative_14[k]), 4))
        avg_heuristic_risk_count_negative_21.append(round (mean(heuristic_risk_count_negative_21[k]), 4))

        avg_heuristic_risk_negative_arrival.append(round (mean(heuristic_risk_negative_arrival[k]), 5))
        avg_heuristic_risk_negative_arrival_7.append(round (mean(heuristic_risk_negative_arrival_7[k]), 5))
        avg_heuristic_risk_negative_arrival_14.append(round (mean(heuristic_risk_negative_arrival_14[k]), 5))
        avg_heuristic_risk_negative_arrival_21.append(round (mean(heuristic_risk_negative_arrival_21[k]), 5))

        total_min_negative_stocks_risk.append(min_negative_stocks_risk[k])
        total_min_negative_stocks_risk_7.append(min_negative_stocks_risk_7[k])
        total_min_negative_stocks_risk_14.append(min_negative_stocks_risk_14[k])
        total_min_negative_stocks_risk_21.append(min_negative_stocks_risk_21[k])

    df.to_csv("result_v2.csv", index=0)

    df = pd.DataFrame(capacity_list, columns= ['Capacity'])

    df['LP model'] = avg_optimal_value
    
    df['Conservative heuristic'] = avg_heuristic_value
    df['Heuristic 1-days-ahead risk'] = avg_heuristic_risk_value
    df['Heuristic 7-days-ahead risk'] = avg_heuristic_risk_value_7
    df['Heuristic 14-days-ahead risk'] = avg_heuristic_risk_value_14
    df['Heuristic 21-days-ahead risk'] = avg_heuristic_risk_value_21
    
    df['LP model stocks'] = avg_stocks_optimal
    
    df['Conservative heuristic stocks'] = avg_stocks_heuristic
    df['Heuristic q-1 risk stocks'] = avg_stocks_risk
    df['Heuristic q-7 risk stocks'] = avg_stocks_risk_7
    df['Heuristic q-14 risk stocks'] = avg_stocks_risk_14
    df['Heuristic q-21 risk stocks'] = avg_stocks_risk_21

    df['Heuristic 1-days-ahead negative on arrival'] = avg_heuristic_risk_negative_arrival
    df['Heuristic 7-days-ahead negative on arrival'] = avg_heuristic_risk_negative_arrival_7
    df['Heuristic 14-days-ahead negative on arrival'] = avg_heuristic_risk_negative_arrival_14
    df['Heuristic 21-days-ahead negative on arrival'] = avg_heuristic_risk_negative_arrival_21

    df['Heuristic 1-days-ahead negative days'] = avg_heuristic_risk_count_negative
    df['Heuristic 7-days-ahead negative days'] = avg_heuristic_risk_count_negative_7
    df['Heuristic 14-days-ahead negative days'] = avg_heuristic_risk_count_negative_14
    df['Heuristic 21-days-ahead negative days'] = avg_heuristic_risk_count_negative_21
    
    df['LP model utilization'] = avg_utilization_T
    df['Heuristic utilization'] = avg_utilization_T_heuristic
    df['Heuristic 1-days-ahead utilization'] = avg_utilization_T_heuristic_risk
    df['Heuristic 7-days-ahead utilization'] = avg_utilization_T_heuristic_risk_7
    df['Heuristic 14-days-ahead utilization'] = avg_utilization_T_heuristic_risk_14
    df['Heuristic 21-days-ahead utilization'] = avg_utilization_T_heuristic_risk_21

    df['Heuristic 1-days-ahead max negative stocks'] = total_min_negative_stocks_risk
    df['Heuristic 7-days-ahead max negative stocks'] = total_min_negative_stocks_risk_7
    df['Heuristic 14-days-ahead max negative stocks'] = total_min_negative_stocks_risk_14
    df['Heuristic 21-days-ahead max negative stocks'] = total_min_negative_stocks_risk_21


    if LAST_DAY_VACCINES:
        df.to_csv("result_summary_v2_zero_stocks_123.csv", index = 0)
    else:
        df.to_csv("result_summary_v2.csv", index = 0)