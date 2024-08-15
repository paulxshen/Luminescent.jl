from luminescent import finetune, apply_design, load_component, load_solution
# sol = finetune(2)
sol = load_solution(study="inverse_design")
c = load_component(study="inverse_design")
c = apply_design(c, sol)
