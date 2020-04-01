# Feather

Feather: Lightweight Multi-party Updatable Delegated Private Set Intersection
With the growth of cloud computing, the need arises for Private Set Intersection (PSI) protocols that can operate on outsourced data and delegate computation to cloud servers. 
One limitation of existing delegated PSI protocols is that they are all designed for static data and do not allow efficient update on outsourced data. 
Another limitation is that they cannot efficiently support PSI among multiple clients, which is often needed in practice. 
This work presents “Feather”, the first delegated PSI protocol that supports efficient data updates and multi-party PSI computation on outsourced datasets. 
The clients can independently prepare and upload their private data to the cloud once, then delegate the computation an unlimited number of times. The update operation has O(1) communication and computation complexity, and this is achieved without sacrificing PSI efficiency and security. Feather does not use public key cryptography, that makes it more scalable. 
We have implemented a prototype and compared the concrete performance against the state of the art. 
The evaluation indicates that Feather does achieve better performance in both update and PSI computation.
