# def clean_file_parsing(file_object):
# def read_init(file_path):
# def read_solver(file_path):
# def write_init(file_obj, init):
# def write_solver(file_obj, solver):
# def get_cost(estimator_binary, model, solver, estimator, param=None):
# def simulate(simulator, model, init, solver, param=None, init_out=None, 
# def sim_result_to_array(sim_std_out, output_mode=2):
# def read_parameter(parameter_file_location):
# def param_from_opt(opt_file):
# def optimise(optimiser_path, model, solver, estimator, optimiser, *args):

clean_file_parsing = function(file.object) {
  file.object = sapply(file.object, function(string) gsub(pattern = "#[^\\\n]*", replacement = "", x = string))
  file.object = file.object[which(file.object != "")]
  file.object
}

read.init = function(file.path) {
  data = read.table(file.path, skip = 1, blank.lines.skip = TRUE, comment.char = "#")
  colnames(data)[1:4] = c("x","y","z","r")
  data
}

read.solver = function(file.path) {
  # """ Read an organism solver file and return a list of strings 
  # 
  #   :file_path: string specifying where to read the solver file from.
  #   :returns: solver -- a list of strings from the solver file.
  #   """
  # solver=[]
  # with open(file_path, 'r') as solver_file:
  #   for line in clean_file_parsing(solver_file):
  #   solver.append([entry for entry in line.split()])
  # 
  # solver = [ entry for sublist in solver for entry in sublist] # Flatten
  # return solver 
  file.path = "~/projects/sum16_solver/MSB_WCK/solvers/rk4"
  
  
}