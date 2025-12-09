import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#add the parent directory to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.poisson import solve  # Importez votre fonction solve

def validate_poisson_solver():
    """
    Validation complète de votre solveur Poisson
    """
    
    # ==================== TEST 1: Solution analytique connue ====================
    print("="*60)
    print("TEST 1: Validation avec solution analytique sinusoïdale")
    print("="*60)
    
    # Solution analytique: u(x,y) = sin(πx)sin(πy)
    # ∇²u = -2π² sin(πx)sin(πy)
    def f_test1(x, y):
        return -2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # Paramètres
    a, b = 1.0, 1.0  # Domaine [0,1] x [0,1]
    h, k = 0.1, 0.1  # Pas spatial
    n, m = int(1 + a/h), int(1 + b/k)
    
    # Conditions aux limites (u=0 sur tout le bord pour cette solution)
    conditions_lim = [0] * (2*(n+m)-4)
    
    # Calcul solution numérique
    u_numerical = solve(f_test1, a, b, h, k, conditions_lim)
    
    # Calcul solution analytique sur la grille
    x_vals = np.linspace(0, a, n)
    y_vals = np.linspace(0, b, m)
    X, Y = np.meshgrid(x_vals, y_vals)
    u_analytic = np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # Calcul erreurs
    error = u_numerical - u_analytic
    l2_error = np.sqrt(np.mean(error**2))
    max_error = np.max(np.abs(error))
    
    print(f"Grille: {n}x{m} points")
    print(f"Erreur L2: {l2_error:.2e}")
    print(f"Erreur Max: {max_error:.2e}")
    print(f"Valeurs max solution: Numérique={np.max(u_numerical):.4f}, "
          f"Analytique={np.max(u_analytic):.4f}")
    
    # ==================== TEST 2: Solution polynomiale ====================
    print("\n" + "="*60)
    print("TEST 2: Validation avec solution polynomiale")
    print("="*60)
    
    # Solution: u(x,y) = x(1-x)y(1-y)
    # ∇²u = -2[y(1-y) + x(1-x)]
    def f_test2(x, y):
        return -2 * (y*(1-y) + x*(1-x))
    
    # Conditions aux limites (u=0 sur le bord)
    u_numerical2 = solve(f_test2, a, b, h, k, conditions_lim)
    
    # Solution analytique
    u_analytic2 = X * (1-X) * Y * (1-Y)
    error2 = u_numerical2 - u_analytic2
    
    print(f"Erreur L2: {np.sqrt(np.mean(error2**2)):.2e}")
    print(f"Erreur Max: {np.max(np.abs(error2)):.2e}")
    
    # ==================== TEST 3: Vérification symétrie ====================
    print("\n" + "="*60)
    print("TEST 3: Vérification symétrie (f symétrique)")
    print("="*60)
    
    def f_symmetric(x, y):
        return np.exp(-((x-0.5)**2 + (y-0.5)**2))
    
    u_numerical3 = solve(f_symmetric, a, b, h, k, conditions_lim)
    
    # Vérifier symétrie par rapport à x=0.5 et y=0.5
    mid_i, mid_j = m//2, n//2
    diff_symmetry = np.max(np.abs(u_numerical3 - u_numerical3[::-1, :]))  # Symétrie verticale
    diff_symmetry2 = np.max(np.abs(u_numerical3 - u_numerical3[:, ::-1]))  # Symétrie horizontale
    
    print(f"Différence symétrie verticale: {diff_symmetry:.2e}")
    print(f"Différence symétrie horizontale: {diff_symmetry2:.2e}")
    
    # ==================== TEST 4: Test de convergence ====================
    print("\n" + "="*60)
    print("TEST 4: Test d'ordre de convergence")
    print("="*60)
    
    resolutions = [0.2, 0.1, 0.05, 0.025]  # Différents pas
    errors = []
    
    for h_val in resolutions:
        n_val = int(1 + a/h_val)
        m_val = int(1 + b/h_val)
        conditions = [0] * (2*(n_val+m_val)-4)
        
        u_num = solve(f_test1, a, b, h_val, h_val, conditions)
        
        # Grille correspondante pour solution analytique
        x = np.linspace(0, a, n_val)
        y = np.linspace(0, b, m_val)
        Xg, Yg = np.meshgrid(x, y)
        u_ana = np.sin(np.pi * Xg) * np.sin(np.pi * Yg)
        
        errors.append(np.sqrt(np.mean((u_num - u_ana)**2)))
    
    # Calcul taux convergence
    print("Erreurs L2 pour différentes résolutions:")
    for i, (h_val, err) in enumerate(zip(resolutions, errors)):
        print(f"  h={h_val:.3f} -> erreur={err:.2e}")
    
    # Ordre de convergence (devrait être ~2)
    rates = []
    for i in range(1, len(errors)):
        rate = np.log(errors[i-1]/errors[i]) / np.log(2)
        rates.append(rate)
    
    print(f"\nTaux de convergence estimés: {rates}")
    print("(Devrait être proche de 2 pour différences finies d'ordre 2)")
    
    # ==================== VISUALISATION ====================
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Test 1
    im1 = axes[0,0].contourf(X, Y, u_numerical, levels=50, cmap='viridis')
    plt.colorbar(im1, ax=axes[0,0])
    axes[0,0].set_title('Test 1: Solution numérique')
    
    axes[0,1].contourf(X, Y, u_analytic, levels=50, cmap='viridis')
    axes[0,1].set_title('Test 1: Solution analytique')
    
    im3 = axes[0,2].contourf(X, Y, np.abs(error), levels=50, cmap='hot')
    plt.colorbar(im3, ax=axes[0,2])
    axes[0,2].set_title(f'Erreur absolue (max={max_error:.1e})')
    
    # Test 2
    axes[1,0].contourf(X, Y, u_numerical2, levels=50, cmap='viridis')
    axes[1,0].set_title('Test 2: Solution polynomiale')
    
    # Test 3
    axes[1,1].contourf(X, Y, u_numerical3, levels=50, cmap='viridis')
    axes[1,1].set_title('Test 3: Source gaussienne')
    
    # Courbe de convergence
    axes[1,2].loglog(resolutions, errors, 'o-', linewidth=2, markersize=8)
    axes[1,2].loglog(resolutions, [err*4 for err in errors], 'r--', 
                     label='Pente 2 (référence)')
    axes[1,2].set_xlabel('Pas spatial h')
    axes[1,2].set_ylabel('Erreur L2')
    axes[1,2].set_title('Courbe de convergence')
    axes[1,2].legend()
    axes[1,2].grid(True, which="both", ls="--")
    
    plt.tight_layout()
    plt.show()
    
    # ==================== CRITÈRES DE SUCCÈS ====================
    print("\n" + "="*60)
    print("CRITÈRES DE VALIDATION")
    print("="*60)
    
    success_criteria = [
        (l2_error < 1e-3, "Erreur L2 < 1e-3"),
        (max_error < 5e-3, "Erreur Max < 5e-3"),
        (all(r > 1.5 for r in rates), "Ordre convergence > 1.5"),
        (diff_symmetry < 1e-10, "Symétrie préservée")
    ]
    
    all_passed = True
    for condition, description in success_criteria:
        status = "✓ PASS" if condition else "✗ FAIL"
        print(f"{status}: {description}")
        if not condition:
            all_passed = False
    
    if all_passed:
        print("\n✅ VOTRE SOLVEUR POISSON EST VALIDÉ !")
    else:
        print("\n⚠️  Certains tests ont échoué. Vérifiez votre implémentation.")
    
    return u_numerical, u_analytic, error

if __name__ == "__main__":
    # Exécuter la validation
    u_num, u_ana, err = validate_poisson_solver()
    
    # Test additionnel: Vérifier conservation avec source nulle
    print("\n" + "="*60)
    print("TEST ADDITIONNEL: Source nulle (devrait donner solution nulle)")
    print("="*60)
    
    def f_zero(x, y):
        return 0
    
    a, b, h = 1.0, 1.0, 0.1
    n, m = int(1+a/h), int(1+b/h)
    conditions = [1] * (2*(n+m)-4)  # Conditions non nulles
    
    u_zero = solve(f_zero, a, b, h, h, conditions)
    
    # Calculer la partie harmonique (devrait être solution de Laplace)
    laplacian_approx = np.zeros_like(u_zero)
    for i in range(1, m-1):
        for j in range(1, n-1):
            laplacian_approx[i,j] = (u_zero[i+1,j] + u_zero[i-1,j] - 2*u_zero[i,j])/h**2 \
                                  + (u_zero[i,j+1] + u_zero[i,j-1] - 2*u_zero[i,j])/h**2
    
    print(f"Max Laplacien (devrait être ~0): {np.max(np.abs(laplacian_approx)):.2e}")