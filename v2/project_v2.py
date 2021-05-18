import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os
from pathlib import Path
from statistics import mean

PRINT = False
DELTA = {'Pfizer': 21, 'Moderna': 28, 'Astrazeneca': 78}
CSV_INPUT_FOLDER = "input_csv"
CSV_OUTPUT_FOLDER = "csv_solution_v2"
T = 180
NUMBER_OF_ELEMENT = 30
LIMIT_CAPACITY = 8
INCREMENT = 0.2

def optimize_test_capacity_multiple_vaccines(T, B, Delta, Capacity):

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
    m.addConstrs( (first_doses[j,i] + first_doses[j,i-Delta[j]] + stocks[j,i] == B[j,i] + stocks[j, i-1] for j, i in combinations if i >= 1 and i >= Delta[j]))

    # m.addConstrs( (stocks[j,i] <= 1000000  for j, i in combinations if i == T-1 ) )

    m.setObjective((second_doses.prod(time_frame) + 1000 * ( stocks['Pfizer', T-1] + stocks['Moderna', T-1] + stocks['Astrazeneca', T-1])), GRB.MINIMIZE)

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

def heuristic_v2_sum(first_doses_list, delta, t, l):
    value_sum = 0

    for i in range(1, l+1):
        if t + i - delta >= 0 and t + i - delta <= 180:
            value_sum += first_doses_list[t + i - delta]
    return value_sum

def heuristic_v2_risk(b_list, capacity, l):

    vaccines = { "Pfizer": [21, [0]*180, [0]*180], "Moderna": [28, [0]*180, [0]*180], "Astrazeneca": [78, [0]*180, [0]*180] }
    stocks = []

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
                    first_doses_somministrated = stocks[t-1] - first_doses[t-delta] - heuristic_v2_sum(first_doses, delta, t, l) + arrival[t]

                if first_doses_somministrated < 0:
                    first_doses_somministrated = 0

                if t-delta >= 0:
                    second_doses_somministrated = first_doses[t-delta]
                else:
                    second_doses_somministrated = 0

                if first_doses_somministrated == 0 and second_doses_somministrated == 0:
                    fraction_of_capacity = 0
                else:
                    if (first_doses_somministrated + second_doses_somministrated) / remain_capacity > 1:
                        fraction_of_capacity = remain_capacity / (first_doses_somministrated + second_doses_somministrated)
                    else:
                        fraction_of_capacity = (first_doses_somministrated + second_doses_somministrated) / remain_capacity

                if t-delta >= 0:
                    second_doses_somministrated = first_doses[t-delta]
                else:
                    second_doses_somministrated = 0

                '''if index == 0:
                    if (fraction_of_capacity > 0.33):
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.33)
                    else:
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * fraction_of_capacity)
                elif index == 1:
                    if (fraction_of_capacity > 0.5):
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.5)
                    else:
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * fraction_of_capacity)
                elif index == 2:
                    first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, remain_capacity)'''

                first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.33)

                total_second_doses += first_doses_somministrated
                object_function += (t+1+delta)*first_doses_somministrated
                vaccines[u][1][t] = first_doses_somministrated
                remain_capacity -= first_doses_somministrated
                
                if t == 0:
                    vaccines[u][2][t] = - first_doses_somministrated  + arrival[t]
                elif t - delta < 0:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated + arrival[t]
                else:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated - first_doses[t-delta] + arrival[t]
            else:
                vaccines[u][2][t] = stocks[t-1]  - first_doses[t-delta]  + arrival[t]

            index += 1

    sum_penality += vaccines["Pfizer"][2][180-1] 
    sum_penality += vaccines["Astrazeneca"][2][180-1] 
    sum_penality += vaccines["Moderna"][2][180-1]

    return ([object_function/total_second_doses, sum_penality])

def heuristic_v2(b_list, capacity):

    vaccines = { "Pfizer": [21, [0]*180, [0]*180], "Moderna": [28, [0]*180, [0]*180], "Astrazeneca": [78, [0]*180, [0]*180] }
    stocks = []

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

                if first_doses_somministrated == 0 and second_doses_somministrated == 0:
                    fraction_of_capacity = 0
                else:
                    if (first_doses_somministrated + second_doses_somministrated) / remain_capacity > 1:
                        fraction_of_capacity = remain_capacity / (first_doses_somministrated + second_doses_somministrated)
                    else:
                        fraction_of_capacity = (first_doses_somministrated + second_doses_somministrated) / remain_capacity


                if t-delta >= 0:
                    second_doses_somministrated = first_doses[t-delta]
                else:
                    second_doses_somministrated = 0
                
                '''if index == 0:
                    if (fraction_of_capacity > 0.33):
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.33)
                    else:
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * fraction_of_capacity)
                elif index == 1:
                    if (fraction_of_capacity > 0.5):
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.5)
                    else:
                        first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * fraction_of_capacity)
                elif index == 2:
                    first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, remain_capacity)'''

                first_doses_somministrated = min(first_doses_somministrated + second_doses_somministrated, capacity * 0.33)

                total_second_doses += first_doses_somministrated
                object_function += (t+1+delta)*first_doses_somministrated
                vaccines[u][1][t] = first_doses_somministrated
                remain_capacity -= first_doses_somministrated                

                if t == 0:
                    vaccines[u][2][t] = - first_doses_somministrated - first_doses_somministrated + arrival[t]
                else:
                    vaccines[u][2][t] = stocks[t-1] - first_doses_somministrated - first_doses_somministrated  + arrival[t]

            else:
                vaccines[u][2][t] = stocks[t-1] + arrival[t]

            index += 1

    sum_penality += vaccines["Pfizer"][2][180-1] 
    sum_penality += vaccines["Astrazeneca"][2][180-1] 
    sum_penality += vaccines["Moderna"][2][180-1]

    return ([object_function/total_second_doses, sum_penality])

if __name__ == "__main__":
    
    penality_optimal_result = {}
    optimal_result = {}
    
    heuristic_result = {}
    heuristic_risk_result = {}
    heuristic_risk_result_7 = {}
    heuristic_risk_result_14 = {}
    heuristic_risk_result_21 = {}

    penality_heuristic = {}
    penality_heuristic_risk = {}
    penality_heuristic_risk_7 = {}
    penality_heuristic_risk_14 = {}
    penality_heuristic_risk_21 = {}

    num = 1
    while(num < LIMIT_CAPACITY):
        penality_optimal_result[str(num) + " c"] = []
        optimal_result[str(num) + " c"] = []
        
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

        num += INCREMENT
        num = round(num, 2)

    file_list = os.listdir(CSV_INPUT_FOLDER)
    instances = np.arange(1, T + 1, 1).tolist()

    for i in range(0, len(os.listdir(CSV_INPUT_FOLDER)) -2  - (1000 - NUMBER_OF_ELEMENT) ):

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

            for j in result[0]:
                df["first_doses_" + j] = result[0][j]
                df["second_doses_" + j] = result[1][j]
                df["stock_values_" + j] = result[2][j]
                second_doses_sum += sum(result[1][j])
                penality_sum += result[2][j][180-1]

            # Build the capacity list to append at the end of the CSV file
            df["capacity"] = [capacity[u]] * 180
            
            # Calculate the optimal result and append to the CSV file
            solution_penality = 1000 * penality_sum
            optimal_result_without_penality = (result[3] - solution_penality) / second_doses_sum
            optimal_result[u].append( optimal_result_without_penality )
            penality_optimal_result[u].append( penality_sum )

            # Calculate the heuristic result
            heuristic_result[u].append( heu_result[0] )
            penality_heuristic[u].append( heu_result[1] )

            heuristic_risk_result[u].append( heu_result_risk[0] )
            penality_heuristic_risk[u].append( heu_result_risk[1] )

            heuristic_risk_result_7[u].append( heu_result_risk_7[0] )
            penality_heuristic_risk_7[u].append( heu_result_risk_7[1] )

            heuristic_risk_result_14[u].append( heu_result_risk_14[0] )
            penality_heuristic_risk_14[u].append( heu_result_risk_14[1] )

            heuristic_risk_result_21[u].append( heu_result_risk_21[0] )
            penality_heuristic_risk_21[u].append( heu_result_risk_21[1] )

            # Write the solution to the single file (create direcitory if not exists)
            output_dir = Path("csv_solution_v2/capacity_" + u)
            output_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(CSV_OUTPUT_FOLDER + "/capacity_" + u + "/solution_" + file_list[i] )

    instances = np.arange(1, len(optimal_result["1 c"]) + 1 - (1000 - NUMBER_OF_ELEMENT), 1).tolist()

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

    for k in optimal_result:
        df[k + " - optimal solution"] = optimal_result[k]
        df[k + " - optimal remain stocks"] = penality_optimal_result[k]

        df[k + " - heuristic solution"] = heuristic_result[k]
        df[k + " - heuristic remain stocks"] = penality_heuristic[k]

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

        avg_result_difference.append(round( (mean(heuristic_result[k]) - mean(optimal_result[k])) / mean(heuristic_result[k]), 4))
        avg_result_difference_risk.append(round( (mean(heuristic_risk_result[k]) - mean(optimal_result[k])) / mean(heuristic_risk_result[k]), 4))
        avg_result_difference_risk_7.append(round( (mean(heuristic_risk_result_7[k]) - mean(optimal_result[k])) / mean(heuristic_risk_result_7[k]), 4))
        avg_result_difference_risk_14.append(round( (mean(heuristic_risk_result_14[k]) - mean(optimal_result[k])) / mean(heuristic_risk_result_14[k]), 4))
        avg_result_difference_risk_21.append(round( (mean(heuristic_risk_result_21[k]) - mean(optimal_result[k])) / mean(heuristic_risk_result_21[k]), 4))

        avg_stocks_optimal.append(round( mean(penality_optimal_result[k] ), 2))
        avg_stocks_heuristic.append(round( mean(penality_heuristic[k] ), 2))
        avg_stocks_risk.append(round( mean(penality_heuristic_risk[k] ), 2))
        avg_stocks_risk_7.append(round( mean(penality_heuristic_risk_7[k] ), 2))
        avg_stocks_risk_14.append(round( mean(penality_heuristic_risk_14[k] ), 2))
        avg_stocks_risk_21.append(round( mean(penality_heuristic_risk_21[k] ), 2))

        avg_stocks_difference.append(round( (mean(penality_heuristic[k]) - mean(penality_optimal_result[k])) / mean(penality_heuristic[k]), 4))
        avg_stocks_difference_risk.append(round( (mean(penality_heuristic_risk[k]) - mean(penality_optimal_result[k])) / mean(penality_heuristic_risk[k]), 4))
        avg_stocks_difference_risk_7.append(round( (mean(penality_heuristic_risk_7[k]) - mean(penality_optimal_result[k])) / mean(penality_heuristic_risk_7[k]), 4))
        avg_stocks_difference_risk_14.append(round( (mean(penality_heuristic_risk_14[k]) - mean(penality_optimal_result[k])) / mean(penality_heuristic_risk_14[k]), 4))
        avg_stocks_difference_risk_21.append(round( (mean(penality_heuristic_risk_21[k]) - mean(penality_optimal_result[k])) / mean(penality_heuristic_risk_21[k]), 4))

    df.to_csv("result_v2.csv", index=0)

    df = pd.DataFrame(capacity_list, columns= ['Capacity'])

    df['Optimal value'] = avg_optimal_value
    
    df['Heuristic value'] = avg_heuristic_value
    df['Heuristic risk value'] = avg_heuristic_risk_value
    df['Heuristic 7 risk value'] = avg_heuristic_risk_value_7
    df['Heuristic 14 risk value'] = avg_heuristic_risk_value_14
    df['Heuristic 21 risk value'] = avg_heuristic_risk_value_21
    
    df['Result difference'] = avg_result_difference
    df['Heuristic q-1 risk difference'] = avg_result_difference_risk
    df['Heuristic q-7 risk difference'] = avg_result_difference_risk_7
    df['Heuristic q-14 risk difference'] = avg_result_difference_risk_14
    df['Heuristic q-21 risk difference'] = avg_result_difference_risk_21

    df['Optimal rem. stocks'] = avg_stocks_optimal
    
    df['Heuristic rem. stocks'] = avg_stocks_heuristic
    df['Heuristic q-1 risk rem.'] = avg_stocks_risk
    df['Heuristic q-7 risk rem.'] = avg_stocks_risk_7
    df['Heuristic q-14 risk rem.'] = avg_stocks_risk_14
    df['Heuristic q-21 ris rem. stocks'] = avg_stocks_risk_21

    df['Stocks difference'] = avg_stocks_difference
    df['Stocks q-1 risk difference'] = avg_stocks_difference_risk
    df['Stocks q-7 risk difference'] = avg_stocks_difference_risk_7
    df['Stocks q-14 risk difference'] = avg_stocks_difference_risk_14
    df['Stocks q-21 risk difference'] = avg_stocks_difference_risk_21

    df.to_csv("result_summary_v2.csv", index = 0)
