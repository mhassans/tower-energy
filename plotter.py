import sys
import yaml

def towersPerStage1():
    

def main():
    try:
        config_file = sys.argv[1]
    except IndexError:
        print("Please give a valid config file")
        exit()
    try:
        with open(config_file, 'r') as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
    except EnvironmentError:
        print("Please give a valid config file")
        exit()
    
    if (config['function']['towersPerStage1']):
        towersPerStage1()

if __name__ == "__main__":
    start = time.time()
    main()
    print('The program ran in', time.time() - start, 'seconds!')
