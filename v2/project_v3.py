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
NUMBER_OF_ELEMENT = 5
LIMIT_CAPACITY = 1
INCREMENT = 0.05
LAST_DAY_VACCINES = True
PERCENTAGE = 0

def optimize_test_capacity_multiple_vaccines(T, B, real_B, time_delta, Delta, Capacity, start_time, remain_stocks, planned_first_doses):

    m = gp.Model("vaccinations")

    m.Params.LogToConsole = 0

    dicti = {}
    dict_B = {}
    dict_real_B = {}

    for i in B:
        for j in range(0, T):
            dicti[(i, j)] = [j+1]
            dict_B[(i,j)] = B[i][j]  
            dict_real_B[(i,j)] = real_B[i][j]
    
    original_B = B
    original_real_B = real_B

    B = dict_B
    real_B = dict_real_B

    combinations, time_frame = gp.multidict(dicti)

    first_doses = m.addVars(combinations, lb=0.0, vtype=GRB.INTEGER, name="First_Doses")
    second_doses = m.addVars(combinations, lb = 0.0, vtype=GRB.INTEGER, name="Second_Doses")
    stocks = m.addVars(combinations, lb=0.0, vtype=GRB.INTEGER, name="Stocks")
    
    m.addConstrs( (first_doses.sum('*',j) + second_doses.sum('*', j)) <= Capacity for j in range(start_time, T))

    m.addConstrs( (first_doses[j, i] == second_doses[j,i+Delta[j]] for j, i in combinations if i < T-Delta[j] ))

    m.addConstrs( (first_doses[j, i] == 0 for j, i in combinations if i >= T-Delta[j] ))

    m.addConstrs( (first_doses[j, i] == planned_first_doses[j][i] for j, i in combinations if i < start_time ))

    m.addConstrs( second_doses[j,i] == 0 for j, i in combinations if i < Delta[j])

    # real_B list
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i] == real_B[j,i] + 0 for j, i in combinations if i == 0 and i < start_time + time_delta))
        # Not start time
    m.addConstrs( (first_doses[j,i] + second_doses[j,i] + stocks[j,i] == real_B[j,i] + stocks[j,i-1] for j, i in combinations if i > 0 and i != start_time and i < start_time + time_delta))    
        # Start time
    m.addConstrs( (first_doses[j,i] + second_doses[j,i] + stocks[j,i] == real_B[j,i] + remain_stocks[j][i-1] for j, i in combinations if i > 0 and i == start_time and i < start_time + time_delta))

    # B list
    m.addConstrs( (first_doses[j,i] + 0 + stocks[j,i] == B[j,i] + 0 for j, i in combinations if i == 0 and i >= start_time + time_delta))
        # Not start time
    m.addConstrs( (first_doses[j,i] + second_doses[j,i] + stocks[j,i] == B[j,i] + stocks[j,i-1] for j, i in combinations if i > 0 and i != start_time and i >= start_time + time_delta))
        # Start time
    m.addConstrs( (first_doses[j,i] + second_doses[j,i] + stocks[j,i] == B[j,i] + remain_stocks[j][i-1] for j, i in combinations if i > 0 and i == start_time and i >= start_time + time_delta))

    m.setObjective(  (gp.quicksum(first_doses[j,i] * (i+Delta[j]) for j,i in combinations) + 1000 * ( stocks['Pfizer', T-1] + stocks['Moderna', T-1] + stocks['Astrazeneca', T-1])), GRB.MINIMIZE)

    m.optimize()

    if (m.solCount > 0):

        resultList = m.getAttr(GRB.Attr.X, m.getVars())
            
        first_doses_dict = {}
        second_doses_dict = {}
        stocks_dict = {}

        for i in range(0, len(original_B)):
            
            first_doses_dict[list(original_B)[i]] = resultList[(T)*i:(T)*(i+1)]
            first_doses_values = resultList[(T)*i:(T)*(i+1)]

        for i in range(0, len(original_B)):

            second_doses_dict[list(original_B)[i]] = resultList[(T)*(i+len(original_B)):(T)*(i+1+len(original_B))]
            second_doses_values = resultList[(T)*(i+len(original_B)):(T)*(i+1+len(original_B))]
        
        for i in range(0, len(original_B)):
            stocks_dict[list(original_B)[i]] = resultList[(T)*(i+2*len(original_B)):(T)*(i+1+2*len(original_B))]
            stock_values = resultList[(T)*(i+2*len(original_B)):(T)*(i+1+2*len(original_B))]

        object_function_value = m.objVal
        
        return[first_doses_dict, second_doses_dict, stocks_dict, object_function_value]
    else:
        print("\n***** No solutions found *****")
        return -1

def dict_to_list(dict):
    value_list = []
    for i in dict:
        value_list.append(dict[i])
    return value_list

def dict_to_mean(dict):
    for i in dict:
        dict[i] = mean(dict[i])
    return dict

def get_column_from_df(df, column_name):
    return df[column_name].values.tolist()

def add_column_to_df(df, values_list, new_column_name):
    df[new_column_name] = values_list
    return df

if __name__ == "__main__":

    shutil.rmtree("csv_solution_v2")
    os.mkdir("csv_solution_v2")

    file_list = os.listdir(CSV_INPUT_FOLDER)
    instances = np.arange(1, T + 1, 1).tolist()

    instances_error = []
    result_list = []

    percentage_list = {}
    list_index = []

    num = 0
    while(num < LIMIT_CAPACITY):
        percentage_list[str(num)] = []
        num += INCREMENT
        num = round(num, 2)
        list_index.append(num)
        
    for i in range(0, len(os.listdir(CSV_INPUT_FOLDER)) - 900):

        print("Processing instance: " + str(i))

        prime_dosi_effettuate = {"Astrazeneca":[0]*180, "Pfizer":[0]*180, "Moderna":[0]*180}
        seconde_dosi_effettuate = {"Astrazeneca":[0]*180, "Pfizer":[0]*180, "Moderna":[0]*180}
        prime_dosi_programmate = {"Astrazeneca":[0]*180, "Pfizer":[0]*180, "Moderna":[0]*180}
        scorte = {"Astrazeneca":[0]*180, "Pfizer":[0]*180, "Moderna":[0]*180}

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
        while(num < 2):
            capacity[num] =  int(num * int(total_capacity/180))
            num += 1
            num = round(num, 2)

        # For each capacity...

        capacity_list = []

        for p in percentage_list:

            for u in capacity:

                capacity_list.append(u)

                penality_sum = 0
                second_doses_sum = 0
                new_doses = 0
                REMAIN_STOCKS = {"Pfizer": 0, "Moderna" : 0, "Astrazeneca" : 0}

                # Calcolo la prima soluzione
                result = optimize_test_capacity_multiple_vaccines(len(b_list_0), b_list, b_list, 0, DELTA, capacity[u], 0, scorte, prime_dosi_effettuate)

                # Se ho errore su un'istanza, la aggiungo alla lista per effettuare un controllo successivamente
                if result == -1:
                    instances_error.append(i)
                    continue
                
                scorte = result[2]
                prime_dosi_programmate = result[0]

                total_value = 0

                for t in range (0, T):

                    # Se ho aggiornamenti, ricalcolo la soluzione ottima
                    if t  > 1000 and t % 20 == 0: # Aggiungere il controllo in caso di ricalcolo
                        scorte["Pfizer"][t-1] += REMAIN_STOCKS["Pfizer"]
                        scorte["Moderna"][t-1] += REMAIN_STOCKS["Moderna"]
                        scorte["Astrazeneca"][t-1] += REMAIN_STOCKS["Astrazeneca"]

                        REMAIN_STOCKS["Pfizer"] = 0
                        REMAIN_STOCKS["Moderna"] = 0
                        REMAIN_STOCKS["Astrazeneca"] = 0

                        result = optimize_test_capacity_multiple_vaccines(len(b_list_0), b_list, b_list, 0, DELTA, capacity[u], t, scorte, prime_dosi_effettuate)

                        # Se ho errore su un'istanza, la aggiungo alla lista per effettuare un controllo successivamente
                        if result == -1:
                            instances_error.append(i)
                            input()
                            continue

                        scorte = result[2]
                        prime_dosi_programmate = result[0]

                    # Aggiorno il valore delle prime dosi somministrate e della soluzione
                    for i in prime_dosi_effettuate:
                        
                        prime_dosi_effettuate[i][t] = int( prime_dosi_programmate[i][t] * (1-float(p)) )
                        REMAIN_STOCKS[i] += prime_dosi_programmate[i][t] - prime_dosi_effettuate[i][t]

                        if (t + DELTA[i] < T):
                            seconde_dosi_effettuate[i][t+DELTA[i]] = prime_dosi_effettuate[i][t]
                            total_value += seconde_dosi_effettuate[i][t+DELTA[i]] * (t+DELTA[i])
                
                # Aggiungo le scorte finali al computo del valore della soluzione
                scorte["Pfizer"][len(scorte["Pfizer"])-1] += REMAIN_STOCKS["Pfizer"]
                scorte["Moderna"][len(scorte["Moderna"])-1] += REMAIN_STOCKS["Moderna"]
                scorte["Astrazeneca"][len(scorte["Astrazeneca"])-1] += REMAIN_STOCKS["Astrazeneca"]      

                second_doses_sum =sum(seconde_dosi_effettuate["Pfizer"]) + sum(seconde_dosi_effettuate["Moderna"]) + sum(seconde_dosi_effettuate["Astrazeneca"])            
                penality_sum = scorte["Pfizer"][len(scorte["Pfizer"])-1] + scorte["Moderna"][len(scorte["Moderna"])-1] + scorte["Astrazeneca"][len(scorte["Astrazeneca"])-1]

                new_doses = 0
                new_doses += (scorte["Pfizer"][len(scorte["Pfizer"])-1] / 2) * (180 + DELTA["Pfizer"])
                new_doses += (scorte["Moderna"][len(scorte["Moderna"])-1] / 2) * (180 + DELTA["Moderna"])
                new_doses += (scorte["Astrazeneca"][len(scorte["Astrazeneca"])-1] / 2) * (180 + DELTA["Astrazeneca"])

                # Add all vaccine last day
                optimal_result_without_penality = (total_value + new_doses) / (second_doses_sum + penality_sum/2)

                # print(optimal_result_without_penality)
                result_list.append(optimal_result_without_penality)

                percentage_list[p].append(optimal_result_without_penality)

    df = pd.DataFrame(list_index, columns= ['Stocks_Percentae'])
    df = df.reset_index()

    percentage_list = dict_to_mean(percentage_list)
    add_column_to_df(df, dict_to_list(percentage_list), "Test")
    df.to_csv("Test.csv", index=0)

    print(instances_error)
    print(mean(result_list))