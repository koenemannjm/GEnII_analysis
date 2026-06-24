import os

def ReadLogFile(logfile):
    with open(logfile,'r') as f:
        lines = f.readlines()

    return lines

def FindCutSummary(logfile):

    statement = "Cut summary:"
    
    lines = ReadLogFile(logfile)
    
    table_lines = []
    find_table_line = 0
    table_index = 0
    for line in lines:
    
        if line.rstrip() == statement:
            table_index = find_table_line

        if table_index > 0 and find_table_line > table_index:
            table_line = line.rstrip()
            #print(table_line)
            table_lines.append(table_line)
            
        find_table_line += 1

    return table_lines

def FindPhysicsEvents(logfile):

    statement = "Counter summary:"
    decoder = "physics events"
    len_decoder = len(decoder)
    
    lines = ReadLogFile(logfile)
    
    physics_events = 0
    find_table_line = 0
    table_index = 0
    for line in lines:
    
        if line.rstrip() == statement:
            table_index = find_table_line

        if table_index > 0 and find_table_line > table_index:
            table_line = line.rstrip()
            #print(table_line)

            for i in range(len(table_line)-len_decoder-1):
                if table_line[i:i+len_decoder] == decoder:

                    physics_events_lines = table_line[:i]
                    physics_events = float(physics_events_lines)
            
        find_table_line += 1

    return physics_events

def GetCutSummaryInfo(logfile, eventinfo, timeinfo):

    columnfinder = "CoarseReconstruct_master"
    len_cfinder = len(columnfinder)
    table_lines = FindCutSummary(logfile)
    physics_events = str(FindPhysicsEvents(logfile))[:-2].rstrip()
    # print(physics_events)
    len_physics_events = len(physics_events)

    location_column = 0
    found_first_column = 0
    for line in table_lines:

        rline = line.rstrip()
        len_rline = len(rline)

        if len_cfinder > len_rline:
            continue
        pot_line = rline[:len_cfinder]
        if pot_line == columnfinder:
            # print(rline)
            for i in range(len(rline[len_cfinder:]) - len_physics_events):
                index = len_cfinder+i
                pot_line2 = rline[index:index+len_physics_events].rstrip()
                # print(pot_line2)
                # print(len(pot_line2))
                if pot_line2 == physics_events and found_first_column == 0:
                    location_column = index
                    found_first_column = 1
                    # print(rline[location_column:])


    results = {}
    for einfo in eventinfo:
        
        len_einfo = len(einfo)
        
        for line in table_lines:
            rline = line.rstrip()
            len_rline = len(rline)
            
            if len_einfo > len_rline:
                continue
            
            if rline[:len_einfo] == einfo:

                # print(rline)

                eline = rline[location_column:].rstrip()
                value_str = eline.split()[1]
                value = float(value_str)
                # print(value)
                results[einfo] = value

    
    for tinfo in timeinfo:
        
        len_tinfo = len(tinfo)

        results[tinfo] = {}

        found_timing_sum = 0
        for line in table_lines:
            rline = line.rstrip()
            len_rline = len(rline)

            if rline == "Timing summary:":
                found_timing_sum = 1

            if found_timing_sum == 1:
                if len_tinfo > len_rline:
                    continue

                if rline[:len_tinfo] == tinfo:
                    # print(rline)
                    tline = rline[len_tinfo:].rstrip()

                    real_time = float(tline.split("Real Time = ")[1].split(" seconds")[0])
                    cpu_time = float(tline.split("Cpu Time = ")[1].split(" seconds")[0])
                    
                    results[tinfo]["rt"] = real_time
                    results[tinfo]["ct"] = cpu_time

    return results

table_event_info = [
    "GoodTrack",
    "GoodBBCAL",
    "GoodSBSTrack",
    "Physics_master",
]

table_time_info = [
    "Trackinig",
    "RawDecode",
    "Decode",
    "Total"
]

file_path = "/work/halla//sbs/sbs-gen/GEN_REPLAYS/pass3/test_SBSGEMS/GEN4/He3/logs/"

print_statement = ""

print_statement += f"File path = {file_path} \n"

list_files = os.listdir(file_path)

col_widths = [12, 10, 12, 12, 14, 16]
col_names = ["run number", "segment", "Tot Events", "SBS Events", "Tot Time (s)", " Rate (event/s)"]
for i in range(len(col_widths)):
    print_statement += f"{col_names[i]:<{col_widths[i]}}"

print_statement += "\n"
print_statement += "-" * sum(col_widths) + "\n"
                             
for file_name in list_files:

    run_number = file_name.split("_")[1]
    segment = file_name.split("_")[4]
    log_path = file_path + file_name

    nevents = FindPhysicsEvents(log_path)
    results = GetCutSummaryInfo(log_path, table_event_info, table_time_info)

    print_statement += f"{run_number:<{col_widths[0]}} {segment:<{col_widths[1]}} {nevents:<{col_widths[2]}.0f} {results['GoodSBSTrack']:<{col_widths[3]}} {results['Total']['rt']:<{col_widths[4]}} {nevents/results['Total']['rt']:<{col_widths[5]}.2f} \n"

print(print_statement)
