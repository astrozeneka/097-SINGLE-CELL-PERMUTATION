"use client";

import { useState } from "react";

interface FileItem {
  id: number;
  file: File | null;
}

export default function NearestCellDistance() {
    
    const [fileItems, setFileItems] = useState<FileItem[]>([
        { id: 1, file: null }
    ]);

    return (<>
    
    </>
    )
}