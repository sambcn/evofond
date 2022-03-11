Hydrogrammes
============

L'objectif de cette partie du programme est d'obtenir le débit liquide en amont au cours d'une période de temps. Pour obtenir cette donnée l'utilisateur du logiciel possède 3 façons de procéder, décrites dans chacune des parties suivantes.

.. toctree::

   ./hydrogramme_lavabre.rst
   ./hydrogramme_interpole.rst
   ./hydrogramme_externe.rst



Hydrogramme Lavabre
===================

Cette méthode permet à l'utilisateur de créer un hydrogramme en renseignant seulement les données clées, à savoir :

* Qmax = le débit de pointe (c'est-à-dire le maximum atteint pendant l'intervalle de temps considéré). Noté *Qmax* dans le code.
* Qbase = débit minimal. Noté *Qbase* dans le code.
* tm = le moment où le débit de pointe est atteint. Noté *tm* dans le code.
* alpha = paramètre de forme (lié au degré des polynomes régissant l'hydrogramme). Nota *alpha* dans le code.
* d = durée de l'hydrogramme. Noté *d* dans le code.
* dt = pas de temps pour la discrétisation. Noté *dt* dans le code.

.. WARNING::

   Dans le code on a en fait un intervalle de temps de durée d+1. À corriger ?

L'histogramme est alors construit avec la formule suivante :

:math:`\forall 0 \leq t \leq d \hspace{1cm} Q(t)=Q_{base}+(Q_{max}-Q_{base})\frac{2(\frac{t}{t_m})^\alpha}{1+(\frac{t}{t_m})^{2\alpha}}`

Ce calcul est fait par la fonction *hydcrue* (à l'aide du module numpy). Cette fonction renvoie doncc la liste des instants discrétisant l'intervalle de temps (*tlav*) ainsi que la liste des débits associés *Q*.

L'appel à *hydcrue* est fait par l'intermédiaire de la fonction *hydrogramme* qui récupère les paramètres listés ci-dessus. Elle effectue aussi un post-traitement sur le couple *tlav, Q* pour obtenir le volume cumulé en intégrant le débit au cours du temps (fonction *integrate* du module scipy).


Hydrogramme par interpolation
=============================

Pour cette méthode, on a un appel intermédiaire à *interpolation* qui récupère des couples (instant, débit) dans les listes t et Q, ainsi que la valeur du pas de temps dans la variable *dt*, fournis par l'utilisateur.

Ensuite, un appel à la fonction *interpol* permet de récupérer l'hydrogramme complété de la façon suivante : on garde tous les points renseignés par l'utilisateur et entre chacun de ces instants on complète l'intervalle de temps par pas de dt. Pour chacun de ces nouveaux abscisses on procède à une interpolation linéaire par morceaux (cf `numpy.interp <https://numpy.org/doc/stable/reference/generated/numpy.interp.html>`)


Hydrogramme externe
===================
